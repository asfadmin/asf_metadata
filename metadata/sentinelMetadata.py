#!/usr/bin/python3

import os
import zipfile
from datetime import datetime, timedelta
import lxml.etree as et
from osgeo import gdal, osr, ogr
import pandas as pd
import geopandas as gpd
import numpy as np
import scipy.interpolate
from shapely import wkt
import requests
import json
import xmltodict


def readZipFile(inFile):

  zf = zipfile.ZipFile(inFile, 'r')
  nameList = zf.namelist()

  return (zf, nameList)


def xml2string(inFile, zf=None):

  if zf == None:
    parser = et.XMLParser(remove_blank_text=True)
    doc = et.parse(inFile, parser)
    string = et.tostring(doc)
  else:
    string = zf.read(inFile)

  return string


def get_latlon_extent(filename):

  src = gdal.Open(filename)
  ulx, xres, _, uly, _, yres  = src.GetGeoTransform()
  lrx = ulx + (src.RasterXSize * xres)
  lry = uly + (src.RasterYSize * yres)

  source = osr.SpatialReference()
  source.ImportFromWkt(src.GetProjection())
  target = osr.SpatialReference()
  target.ImportFromEPSG(4326)
  transform = osr.CoordinateTransformation(source, target)

  (lon1, lat1, _) = transform.TransformPoint(ulx, uly)
  (lon2, lat2, _) = transform.TransformPoint(lrx, uly)
  (lon3, lat3, _) = transform.TransformPoint(ulx, lry)
  (lon4, lat4, _) = transform.TransformPoint(lrx, lry)

  lat_min = min(lat1,lat2,lat3,lat4)
  lat_max = max(lat1,lat2,lat3,lat4)
  lon_min = min(lon1,lon2,lon3,lon4)
  lon_max = max(lon1,lon2,lon3,lon4)

  return lat_min, lat_max, lon_min, lon_max


def extractMetadata(safeDir, zf=None):

  manifest = readManifestFile(os.path.join(safeDir, 'manifest.safe'), zf)
  sentinelFiles = getSentinelFiles(safeDir, manifest)
  for sentinelFile in sentinelFiles:
    if sentinelFile['type'] == 'annotation':
      if sentinelFile['polarization'] == 'VV' or \
        sentinelFile['polarization'] == 'HH':
        if sentinelFile['swath'] == 'IW1':
          annotation_iw1 = readAnnotationFile(sentinelFile['name'], zf)
        elif sentinelFile['swath'] == 'IW2':
          annotation_iw2 = readAnnotationFile(sentinelFile['name'], zf)
        elif sentinelFile['swath'] == 'IW3':
          annotation_iw3 = readAnnotationFile(sentinelFile['name'], zf)

  return (manifest, annotation_iw1, annotation_iw2, annotation_iw3)


def readManifestFile(manifestFile, zf=None):

  ns_xfdu = {'xfdu': 'urn:ccsds:schema:xfdu:1'}
  ns = dict(
    list(ns_xfdu.items())
  )

  ### Read manifest file
  doc = et.fromstring(xml2string(manifestFile, zf))

  ### Pack first level into strings
  meta = {}
  meta['informationPackageMap'] = et.tostring( \
    doc.xpath('/xfdu:XFDU/informationPackageMap', namespaces=ns)[0])
  meta['metadataSection'] = \
    et.tostring(doc.xpath('/xfdu:XFDU/metadataSection', namespaces=ns)[0])
  meta['dataObjectSection'] = \
    et.tostring(doc.xpath('/xfdu:XFDU/dataObjectSection', namespaces=ns)[0])

  return meta


def readAnnotationFile(annotationFile, zf=None):
  
  ### Read annotation file
  doc = et.fromstring(xml2string(annotationFile, zf))

  ### Pack first level into strings
  meta = {}
  meta['adsHeader'] = et.tostring(doc.xpath('/product/adsHeader')[0])
  meta['qualityInformation'] = \
    et.tostring(doc.xpath('/product/qualityInformation')[0])
  meta['generalAnnotation'] = \
    et.tostring(doc.xpath('/product/generalAnnotation')[0])
  meta['imageAnnotation'] = \
    et.tostring(doc.xpath('/product/imageAnnotation')[0])
  meta['dopplerCentroid'] = \
    et.tostring(doc.xpath('/product/dopplerCentroid')[0])
  meta['antennaPattern'] = et.tostring(doc.xpath('/product/antennaPattern')[0])
  meta['swathTiming'] = et.tostring(doc.xpath('/product/swathTiming')[0])
  meta['geolocationGrid'] = \
    et.tostring(doc.xpath('/product/geolocationGrid')[0])
  meta['coordinateConversion'] = \
    et.tostring(doc.xpath('/product/coordinateConversion')[0])
  meta['swathMerging'] = et.tostring(doc.xpath('/product/swathMerging')[0])

  return meta


def annotation2spline(meta):

  ### Extract geolocation grid from metadata
  docGrid = et.fromstring(meta['geolocationGrid'])
  numberOfPoints = int(docGrid.xpath('/geolocationGrid/' \
    'geolocationGridPointList/@count')[0])

  ### Determine line/pixel indices for spline calculation
  lines = []
  pixels = []
  for point in range(numberOfPoints):
    xml = ('/geolocationGrid/geolocationGridPointList/' \
      'geolocationGridPoint[{0}]/line'.format(point+1))
    line = int(docGrid.xpath(xml)[0].text)
    xml = ('/geolocationGrid/geolocationGridPointList/' \
      'geolocationGridPoint[{0}]/pixel'.format(point+1))
    pixel = int(docGrid.xpath(xml)[0].text)
    if line == 0:
      pixels.append(pixel)
    if pixel == 0:
      lines.append(line)
  lineVector = np.array(lines)
  lineCount = len(lines)
  pixelVector = np.array(pixels)
  pixelCount = len(pixels)

  ### Build latitude/longitude grid
  latGrid = np.zeros((lineCount, pixelCount), dtype=float)
  lonGrid = np.zeros((lineCount, pixelCount), dtype=float)
  heightGrid = np.zeros((lineCount, pixelCount), dtype=float)
  for point in range(numberOfPoints):
    xml = ('/geolocationGrid/geolocationGridPointList/' \
      'geolocationGridPoint[{0}]/line'.format(point+1))
    line = int(docGrid.xpath(xml)[0].text)
    xml = ('/geolocationGrid/geolocationGridPointList/' \
      'geolocationGridPoint[{0}]/pixel'.format(point+1))
    pixel = int(docGrid.xpath(xml)[0].text)
    xml = ('/geolocationGrid/geolocationGridPointList/' \
      'geolocationGridPoint[{0}]/latitude'.format(point+1))
    lat = float(docGrid.xpath(xml)[0].text)
    xml = ('/geolocationGrid/geolocationGridPointList/' \
      'geolocationGridPoint[{0}]/longitude'.format(point+1))
    lon = float(docGrid.xpath(xml)[0].text)
    xml = ('/geolocationGrid/geolocationGridPointList/' \
      'geolocationGridPoint[{0}]/height'.format(point+1))
    height = float(docGrid.xpath(xml)[0].text)
    indexLine = lines.index(line)
    indexPixel = pixels.index(pixel)
    latGrid[indexLine, indexPixel] = lat
    lonGrid[indexLine, indexPixel] = lon
    heightGrid[indexLine, indexPixel] = height

  ### Calculate spline for interpolation
  latSpline = scipy.interpolate.RectBivariateSpline(lineVector, pixelVector,
    latGrid)
  lonSpline = scipy.interpolate.RectBivariateSpline(lineVector, pixelVector,
    lonGrid)
  heightSpline = scipy.interpolate.RectBivariateSpline(lineVector, pixelVector,
    heightGrid)

  return (latSpline, lonSpline, heightSpline)


def annotation2burstLocation(meta, burstMap, granuleInfo):

  ### Initialize geodataframe and define data types
  df = pd.DataFrame()

  (latSpline, lonSpline, heightSpline) = annotation2spline(meta)

  ### Extract metadata
  speedOfLight = 299792458.0
  doc = et.fromstring(meta['adsHeader'])
  swath = doc.xpath('/adsHeader/swath')[0].text
  swathNumber = int(swath[-1:])
  polarization = doc.xpath('/adsHeader/polarisation')[0].text
  stopTime = doc.xpath('/adsHeader/stopTime')[0].text

  docGeneral = et.fromstring(meta['generalAnnotation'])
  rangeSamplingRate = float(docGeneral.xpath('/generalAnnotation/' \
    'productInformation/rangeSamplingRate')[0].text)
  radarWavelength = speedOfLight / \
    float(docGeneral.xpath('/generalAnnotation/productInformation/' \
    'radarFrequency')[0].text)
  azimuthSteeringRate = float(docGeneral.xpath('/generalAnnotation/' \
    'productInformation/azimuthSteeringRate')[0].text)
  prf = float(docGeneral.xpath('/generalAnnotation/downlinkInformationList/' \
    'downlinkInformation/prf')[0].text)
  azimuthFmRateListCount = int(docGeneral.xpath('/generalAnnotation/' \
    'azimuthFmRateList/@count')[0])
  fmRates = []
  for rate in range(azimuthFmRateListCount):
    xml = ('/generalAnnotation/azimuthFmRateList/azimuthFmRate[{0}]/' \
      'azimuthTime'.format(rate+1))
    refTime = datetime.fromisoformat(docGeneral.xpath(xml)[0].text)
    xml = ('/generalAnnotation/azimuthFmRateList/azimuthFmRate[{0}]/' \
      'azimuthFmRatePolynomial'.format(rate+1))
    polyCoeffs = list(map(float, docGeneral.xpath(xml)[0].text.split(' ')))
    fmRates.append((refTime, polyCoeffs))

  docImage = et.fromstring(meta['imageAnnotation'])
  azimuthTimeInterval = float(docImage.xpath('/imageAnnotation/' \
    'imageInformation/azimuthTimeInterval')[0].text)
  rangePixelSpacing = float(docImage.xpath('/imageAnnotation/' \
    'imageInformation/rangePixelSpacing')[0].text)
  rangeWindowType = docImage.xpath('/imageAnnotation/processingInformation/' \
    'swathProcParamsList/swathProcParams/rangeProcessing/windowType') \
    [0].text
  rangeWindowCoefficient = float(docImage.xpath('/imageAnnotation/' \
    'processingInformation/swathProcParamsList/swathProcParams/' \
    'rangeProcessing/windowCoefficient')[0].text)
  azimuthWindowType = docImage.xpath('/imageAnnotation/' \
    'processingInformation/swathProcParamsList/swathProcParams/' \
    'azimuthProcessing/windowType')[0].text
  azimuthWindowCoefficient = float(docImage.xpath('/imageAnnotation/' \
    'processingInformation/swathProcParamsList/swathProcParams/' \
    'azimuthProcessing/windowCoefficient')[0].text)

  docDoppler = et.fromstring(meta['dopplerCentroid'])
  dcEstimateListCount = int(docDoppler.xpath('/dopplerCentroid/' \
    'dcEstimateList/@count')[0])
  dopplers = []
  for doppler in range(dcEstimateListCount):
    xml = ('/dopplerCentroid/dcEstimateList/dcEstimate[{0}]/azimuthTime' \
      .format(doppler+1))
    refTime = datetime.fromisoformat(docDoppler.xpath(xml)[0].text)
    xml = ('/dopplerCentroid/dcEstimateList/dcEstimate[{0}]/dataDcPolynomial' \
      .format(doppler+1))
    polyCoeffs = list(map(float, docDoppler.xpath(xml)[0].text.split(' ')))
    dopplers.append((refTime, polyCoeffs))

  docSwath = et.fromstring(meta['swathTiming'])
  linesPerBurst = int(docSwath.xpath('/swathTiming/linesPerBurst')[0].text)
  samplesPerBurst = int(docSwath.xpath('/swathTiming/samplesPerBurst')[0].text)
  numberOfBursts = int(docSwath.xpath('/swathTiming/burstList/@count')[0])

  docGrid = et.fromstring(meta['geolocationGrid'])
  numberOfPoints = int(docGrid.xpath('/geolocationGrid/' \
    'geolocationGridPointList/@count')[0])

  ### Determine line/pixel indices for spline calculation
  lines = []
  pixels = []
  slantRange = []
  for point in range(numberOfPoints):
    xml = ('/geolocationGrid/geolocationGridPointList/' \
      'geolocationGridPoint[{0}]/line'.format(point+1))
    line = int(docGrid.xpath(xml)[0].text)
    xml = ('/geolocationGrid/geolocationGridPointList/' \
      'geolocationGridPoint[{0}]/pixel'.format(point+1))
    pixel = int(docGrid.xpath(xml)[0].text)
    xml = ('/geolocationGrid/geolocationGridPointList/' \
      'geolocationGridPoint[{0}]/slantRangeTime'.format(point+1))
    slantRangeTime = float(docGrid.xpath(xml)[0].text)
    if line == 0:
      pixels.append(pixel)
    if pixel == 0:
      lines.append(line)
      sr = slantRangeTime * speedOfLight / 2.0
      slantRange.append(sr)
  maxPixel = max(pixels)
  centerPixel = maxPixel / 2.0

  acqTime = []
  sensingTime = []
  anxTime = []
  for burst in range(numberOfBursts):
    xml = ('/swathTiming/burstList/burst[{0}]/azimuthAnxTime' \
      .format(burst+1))
    anxTimeBurst = float(docSwath.xpath(xml)[0].text)
    xml = ('/swathTiming/burstList/burst[{0}]/azimuthTime' \
      .format(burst+1))
    azimuthStartTime = docSwath.xpath(xml)[0].text
    xml = ('/swathTiming/burstList/burst[{0}]/sensingTime' \
      .format(burst+1))
    sensingStart = docSwath.xpath(xml)[0].text
    acqTime.append(datetime.fromisoformat(azimuthStartTime))
    sensingTime.append(datetime.fromisoformat(sensingStart))
    anxTime.append(anxTimeBurst)
  acqTime.append(datetime.fromisoformat(stopTime))
  sensingTime.append(datetime.fromisoformat(stopTime))
  sensingPeriod = (linesPerBurst-1)/prf
  sensingMidPeriod = sensingPeriod / 2.0

  ### Determine valid data and calculate burst boundaries
  for burst in range(numberOfBursts):
 
    firstLine = lines[burst]
    lastLine = lines[burst+1]
    startingRange = slantRange[burst]
    burstTimePeriod = (lastLine - firstLine - 1) * azimuthTimeInterval
    acqPeriod = abs(acqTime[burst+1] - acqTime[burst])
    acqTimePeriod = acqPeriod.seconds + acqPeriod.microseconds/1000000
    lastLine += \
      int(abs(burstTimePeriod - acqTimePeriod) / azimuthTimeInterval + 0.5)
    lats = np.zeros(4, dtype=float) 
    lons = np.zeros(4, dtype=float)
    heights = np.zeros(4, dtype=float)
    lats[0] = latSpline(firstLine, 0)[0][0]
    lons[0] = lonSpline(firstLine, 0)[0][0]
    heights[0] = heightSpline(firstLine, 0)[0][0]
    lats[1] = latSpline(firstLine, maxPixel)[0][0]
    lons[1] = lonSpline(firstLine, maxPixel)[0][0]
    heights[1] = heightSpline(firstLine, maxPixel)[0][0]
    lats[2] = latSpline(lastLine, maxPixel)[0][0]
    lons[2] = lonSpline(lastLine, maxPixel)[0][0]
    heights[2] = heightSpline(lastLine, maxPixel)[0][0]
    lats[3] = latSpline(lastLine, 0)[0][0]
    lons[3] = lonSpline(lastLine, 0)[0][0]
    heights[3] = heightSpline(lastLine, 0)[0][0]
    centerLine = firstLine + (lastLine - firstLine) / 2.0
    terrainHeight = heightSpline(centerLine, centerPixel)[0][0]

    xml = ('/swathTiming/burstList/burst[{0}]/azimuthAnxTime' \
      .format(burst+1))
    anxTimeBurst = float(docSwath.xpath(xml)[0].text)
    xml = ('/swathTiming/burstList/burst[{0}]/azimuthTime' \
      .format(burst+1))
    azimuthStartTime = docSwath.xpath(xml)[0].text
    startTime = datetime.fromisoformat(azimuthStartTime)
    stopTime = startTime + timedelta(seconds=burstTimePeriod)
    azimuthStopTime = stopTime.isoformat(timespec='microseconds')
    xml = ('/swathTiming/burstList/burst[{0}]/sensingTime' \
      .format(burst+1))
    sensingStart = docSwath.xpath(xml)[0].text
    startTime = datetime.fromisoformat(sensingStart)
    midTime = startTime + timedelta(seconds=sensingMidPeriod)
    stopTime = startTime + timedelta(seconds=sensingPeriod)
    sensingStop = stopTime.isoformat(timespec='microseconds')

    timeDiff = [np.abs((midTime - val[0]).total_seconds()) for val in fmRates]
    index = np.argmin(timeDiff)
    azimuthFmRatePolynomial = fmRates[index][1]

    timeDiff = [np.abs((midTime - val[0]).total_seconds()) for val in dopplers]
    index = np.argmin(timeDiff)
    dopplerPolynomial = dopplers[index][1]

    xml = ('/swathTiming/burstList/burst[{0}]/firstValidSample' \
      .format(burst+1))
    validSamples = list(map(int, docSwath.xpath(xml)[0].text.split(' ')))
    firstValidSample = np.median(np.array(validSamples)).astype(int)
    validIndices = \
      [i for i, x in enumerate(validSamples) if x == firstValidSample]
    firstValidLine = int(min(validIndices))
    lastValidLine = int(max(validIndices))
    xml = ('/swathTiming/burstList/burst[{0}]/lastValidSample' \
      .format(burst+1))
    validSamples = list(map(int, docSwath.xpath(xml)[0].text.split(' ')))
    lastValidSample = np.median(np.array(validSamples)).astype(int)
    numberValidLines = lastValidLine - firstValidLine - 1
    numberValidSamples = lastValidSample - firstValidSample - 1

    ipf = float(granuleInfo['ipf'])
    if ipf >= 3.4:
      xml = ('/swathTiming/burstList/burst[{0}]/burstId/@absolute' \
        .format(burst+1))
      absoluteESAburstID = int(docSwath.xpath(xml)[0])
      xml = ('/swathTiming/burstList/burst[{0}]/burstId'.format(burst+1))
      relativeESAburstID = int(docSwath.xpath(xml)[0].text) 
    else:
      absoluteESAburstID = -1
      relativeESAburstID = -1

    track = granuleInfo['trackNumber']
    absoluteBurstID = ('%s-%06d' % (getBurstID(burstMap, track, swath, 
      float(anxTime[burst])), granuleInfo['orbitNumber']))
    relativeBurstID = getBurstID(burstMap, track, swath, float(anxTime[burst]))
    if burst > 0:
      previousBurstID = \
        getBurstID(burstMap, track, swath, float(anxTime[burst-1]))
    else:
      previousBurstID = None
    if burst < numberOfBursts - 1:
      nextBurstID = getBurstID(burstMap, track, swath, float(anxTime[burst+1]))
    else:
      nextBurstID = None

    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint_2D(float(lons[0]), float(lats[0]))
    ring.AddPoint_2D(float(lons[1]), float(lats[1]))
    ring.AddPoint_2D(float(lons[2]), float(lats[2]))
    ring.AddPoint_2D(float(lons[3]), float(lats[3]))
    ring.AddPoint_2D(float(lons[0]), float(lats[0]))
    polygon = ogr.Geometry(ogr.wkbPolygon)
    polygon.AddGeometry(ring)

    values = pd.DataFrame()
    values.at[0,'granuleUR'] = ('%s-BURSTS-%s-%s-%02d' \
      % (granuleInfo['granule'], swath, polarization, burst + 1)) 
    values.at[0,'insertDate'] = None
    values.at[0,'productionDateTime'] = datetime.utcnow().isoformat() + 'Z'
    values.at[0,'identifier'] = ('%s-BURSTS' % granuleInfo['granule'])
    values.at[0,'inputGranule1'] = granuleInfo['granule']
    values.at[0,'beginningDateTime'] = azimuthStartTime
    values.at[0,'endingDateTime'] = azimuthStopTime

    values.at[0,'platformShortName'] = granuleInfo['platform']
    values.at[0,'urlOfFrame'] = granuleInfo['urlOfFrame']
    values.at[0,'swath'] = swath
    values.at[0,'startLineInSwath'] = int(firstLine)
    values.at[0,'endLineInSwath'] = int(lastLine)
    values.at[0,'ipf'] = granuleInfo['ipf']

    values.at[0,'numberOfLines'] = int(linesPerBurst)
    values.at[0,'numberOfSamples'] = int(samplesPerBurst)
    values.at[0,'startingRange'] = float(startingRange)
    values.at[0,'sensingStart'] = sensingStart
    values.at[0,'sensingStop'] = sensingStop
    values.at[0,'burstStartUTC'] = azimuthStartTime
    values.at[0,'burstStopUTC'] = azimuthStopTime
    values.at[0,'trackNumber'] = int(granuleInfo['trackNumber'])
    values.at[0,'frameNumber'] = int(-1)
    values.at[0,'orbitNumber'] = int(granuleInfo['orbitNumber'])
    values.at[0,'swathNumber'] = int(swathNumber)
    values.at[0,'burstNumber'] = int(burst) + 1
    values.at[0,'passDirection'] = granuleInfo['passDirection']
    values.at[0,'azimuthSteeringRate'] = float(azimuthSteeringRate)
    values.at[0,'rangePixelSize'] = float(rangePixelSpacing)
    values.at[0,'rangeSamplingRate'] = float(rangeSamplingRate)
    values.at[0,'azimuthTimeInterval'] = float(azimuthTimeInterval)
    values.at[0,'radarWavelength'] = float(radarWavelength)
    values.at[0,'polarization'] = polarization
    values.at[0,'terrainHeight'] = float(terrainHeight)
    values.at[0,'prf'] = float(prf)
    values.at[0,'firstValidLine'] = int(firstValidLine)
    values.at[0,'numValidLines'] = int(numberValidLines)
    values.at[0,'firstValidSample'] = int(firstValidSample)
    values.at[0,'numValidSamples'] = int(numberValidSamples)
    values.at[0,'rangeWindowType'] = rangeWindowType
    values.at[0,'rangeWindowCoefficient'] = float(rangeWindowCoefficient)
    values.at[0,'azimuthWindowType'] = azimuthWindowType
    values.at[0,'azimuthWindowCoefficient'] = float(azimuthWindowCoefficient)
    values.at[0,'azimuthFMRate'] = \
      ','.join(map(str, azimuthFmRatePolynomial))
    values.at[0,'doppler'] = ','.join(map(str, dopplerPolynomial))

    values.at[0,'absoluteBurstID'] = absoluteBurstID
    values.at[0,'relativeBurstID'] = relativeBurstID
    values.at[0,'previousBurstID'] = previousBurstID
    values.at[0,'nextBurstID'] = nextBurstID
    values.at[0,'absoluteESAburstID'] = absoluteESAburstID
    values.at[0,'relativeESAburstID'] = relativeESAburstID
    values.at[0,'timeSinceAnxNode'] = anxTimeBurst
    values.at[0,'ascendingNodeTime'] = granuleInfo['ascendingNodeTime']
    values.at[0,'hasLand'] = True
    values.at[0,'crossesDateLine'] = False
    values.at[0,'landFraction'] = 1.0

    values.at[0,'pointLongitude_1'] = float(lons[0])
    values.at[0,'pointLatitude_1'] = float(lats[0])
    values.at[0,'pointHeight_1'] = float(heights[0])
    values.at[0,'pointLongitude_2'] = float(lons[1])
    values.at[0,'pointLatitude_2'] = float(lats[1])
    values.at[0,'pointHeight_2'] = float(heights[1])
    values.at[0,'pointLongitude_3'] = float(lons[2])
    values.at[0,'pointLatitude_3'] = float(lats[2])
    values.at[0,'pointHeight_3'] = float(heights[2])
    values.at[0,'pointLongitude_4'] = float(lons[3])
    values.at[0,'pointLatitude_4'] = float(lats[3])
    values.at[0,'pointHeight_4'] = float(heights[3])
    values.at[0,'pointLongitude_5'] = float(lons[0])
    values.at[0,'pointLatitude_5'] = float(lats[0])
    values.at[0,'pointHeight_5'] = float(heights[0])
    values.at[0,'geometry'] = wkt.loads(polygon.ExportToWkt())

    df = pd.concat([df, values], ignore_index=True)

  return df


def getBurstID(burstMap, track, swath, anxTime):

  burstTrack = burstMap[burstMap['track'] == track]
  burstSwath = pd.DataFrame(burstTrack[burstTrack['swath'] == swath])
  burstSwath = burstSwath.reset_index()
  deltaAnxTime = 12*24*3600/175
  for ii in range(len(burstSwath)):
    diff = abs(anxTime - burstSwath.at[ii,'anxTime'])
    if diff < deltaAnxTime:
      burstID = burstSwath.at[ii,'burstID']
      deltaAnxTime = diff

  return burstID


def getSentinelBursts(safeDir, burstMapFile, zf=None, urlStr=None):

  ### Read burst map
  print('Reading burst map file ({0}) ...' \
    .format(os.path.basename(burstMapFile)))
  burstMap = gpd.read_file(burstMapFile)

  ### Define namespaces for XML parsing
  ns_safe = {'safe': 'http://www.esa.int/safe/sentinel-1.0'}
  ns_s1 = {'s1': 'http://www.esa.int/safe/sentinel-1.0/sentinel-1'}
  ns = dict(
    list(ns_safe.items()) +
    list(ns_s1.items())
  )

  '''
  (safeDir, manifest, annotation_iw1, annotation_iw2, annotation_iw3) = \
    extractMetadataFromZip(zipFile)
  '''
  (manifest, annotation_iw1, annotation_iw2, annotation_iw3) = \
    extractMetadata(safeDir, zf)
  print('Extracting information for granule ({0}) ...'.format(safeDir[:-6]))

  doc = et.fromstring(manifest['metadataSection'])
  metaObjectCount = len(doc.xpath('/metadataSection/metadataObject/@ID'))
  for ii in range(metaObjectCount):
    id = doc.xpath('/metadataSection/metadataObject[{0}]/@ID' \
      .format(ii+1))[0]
    if id == 'processing':
      ipf = doc.xpath('/metadataSection/metadataObject[{0}]/metadataWrap/' \
        'xmlData/safe:processing/safe:facility/safe:software/@version' \
        .format(ii+1), namespaces=ns)[0]
    elif id == 'platform':
      param = doc.xpath('/metadataSection/metadataObject[{0}]/metadataWrap/' \
        'xmlData/safe:platform/safe:familyName'.format(ii+1), namespaces=ns)
      platform = param[0].text
      param = doc.xpath('/metadataSection/metadataObject[{0}]/metadataWrap/' \
        'xmlData/safe:platform/safe:number'.format(ii+1), namespaces=ns)
      platform += param[0].text
    elif id == 'measurementOrbitReference':
      param = doc.xpath('/metadataSection/metadataObject[{0}]/metadataWrap/' \
        'xmlData/safe:orbitReference/safe:orbitNumber[@type="start"]' \
        .format(ii+1), namespaces=ns)
      absoluteOrbitNumber = int(param[0].text)
      param = doc.xpath('/metadataSection/metadataObject[{0}]/metadataWrap/' \
        'xmlData/safe:orbitReference/safe:relativeOrbitNumber[@type="start"]' \
        .format(ii+1), namespaces=ns)
      relativeOrbitNumber = int(param[0].text)
      param = doc.xpath('/metadataSection/metadataObject[{0}]/metadataWrap/' \
        'xmlData/safe:orbitReference/safe:extension/s1:orbitProperties/' \
        's1:pass'.format(ii+1), namespaces=ns)
      orbitDirection = param[0].text
      param = doc.xpath('/metadataSection/metadataObject[{0}]/metadataWrap/' \
        'xmlData/safe:orbitReference/safe:extension/s1:orbitProperties/' \
        's1:ascendingNodeTime'.format(ii+1), namespaces=ns)
      ascendingNodeTime = param[0].text

  if not urlStr:
    ### Determine URL (if needed)
    url = 'https://api.daac.asf.alaska.edu/services/search/param?granule_list='
    url += safeDir[:-6]
    #print(url)
    data = requests.get(url).content
    data = json.dumps(xmltodict.parse(data.decode()))
    string = json.loads(data)
    urlStr = string['metalink']['files']['file'][1]['resources']['url']['#text']

  ### Split data files into bursts
  granuleInfo = {}
  granuleInfo['platform'] = platform
  granuleInfo['orbitNumber'] = absoluteOrbitNumber
  granuleInfo['relativeOrbitNumber'] = relativeOrbitNumber
  granuleInfo['trackNumber'] = relativeOrbitNumber
  granuleInfo['passDirection'] = orbitDirection
  granuleInfo['ascendingNodeTime'] = ascendingNodeTime
  granuleInfo['granule'] = safeDir[:-6]
  granuleInfo['urlOfFrame'] = urlStr
  granuleInfo['ipf'] = ipf

  #df = initializeGeoDataFrame()
  df = pd.DataFrame()

  ### Extraction burst information for all IW swaths
  burstInfo = annotation2burstLocation(annotation_iw1, burstMap, granuleInfo)
  df = pd.concat([df, burstInfo], ignore_index=True)
  burstInfo = annotation2burstLocation(annotation_iw2, burstMap, granuleInfo)
  df = pd.concat([df, burstInfo], ignore_index=True)
  burstInfo = annotation2burstLocation(annotation_iw3, burstMap, granuleInfo)
  df = pd.concat([df, burstInfo], ignore_index=True)
  df = fixDataTypes(df)
  gdf = gpd.GeoDataFrame(df, crs='EPSG:4326')
  gdf.set_geometry('geometry')

  return gdf


def fixDataTypes(df):

  df['insertDate'] = pd.to_datetime(df['insertDate'])
  df['productionDateTime'] = pd.to_datetime(df['productionDateTime'])
  df['beginningDateTime'] = pd.to_datetime(df['beginningDateTime'])
  df['endingDateTime'] = pd.to_datetime(df['endingDateTime'])

  df['startLineInSwath'] = df['startLineInSwath'].astype('int32')
  df['endLineInSwath'] = df['endLineInSwath'].astype('int32')

  df['numberOfLines'] = df['numberOfLines'].astype('int32')
  df['numberOfSamples'] = df['numberOfSamples'].astype('int32')
  df['startingRange'] = df['startingRange'].astype('float32')
  df['sensingStart'] = pd.to_datetime(df['sensingStart'])
  df['sensingStop'] = pd.to_datetime(df['sensingStop'])
  df['burstStartUTC'] = pd.to_datetime(df['burstStartUTC'])
  df['burstStopUTC'] = pd.to_datetime(df['burstStopUTC'])
  df['trackNumber'] = df['trackNumber'].astype('int32')
  df['frameNumber'] = df['frameNumber'].astype('int32')
  df['orbitNumber'] = df['orbitNumber'].astype('int32')
  df['swathNumber'] = df['swathNumber'].astype('int32')
  df['burstNumber'] = df['burstNumber'].astype('int32')
  df['azimuthSteeringRate'] = df['azimuthSteeringRate'].astype('float32')
  df['rangePixelSize'] = df['rangePixelSize'].astype('float32')
  df['rangeSamplingRate'] = df['rangeSamplingRate'].astype('float32')
  df['azimuthTimeInterval'] = df['azimuthTimeInterval'].astype('float32')
  df['radarWavelength'] = df['radarWavelength'].astype('float32')
  df['terrainHeight'] = df['terrainHeight'].astype('float32')
  df['prf'] = df['prf'].astype('float32')
  df['firstValidLine'] = df['firstValidLine'].astype('int32')
  df['numValidLines'] = df['numValidLines'].astype('int32')
  df['firstValidSample'] = df['firstValidSample'].astype('int32')
  df['numValidSamples'] = df['numValidSamples'].astype('int32')
  df['rangeWindowCoefficient'] = \
    df['rangeWindowCoefficient'].astype('float32')
  df['azimuthWindowCoefficient'] = \
    df['azimuthWindowCoefficient'].astype('float32')
  
  df['absoluteESAburstID'] = df['absoluteESAburstID'].astype('int32')
  df['relativeESAburstID'] = df['relativeESAburstID'].astype('int32')
  df['timeSinceAnxNode'] = df['timeSinceAnxNode'].astype('float32')
  df['ascendingNodeTime'] = pd.to_datetime(df['ascendingNodeTime'])
  df['hasLand'] = df['hasLand'].astype('bool')
  df['crossesDateLine'] = df['crossesDateLine'].astype('bool')
  df['landFraction'] = df['landFraction'].astype('float32')

  df['pointLongitude_1'] = df['pointLongitude_1'].astype('float32')
  df['pointLatitude_1'] = df['pointLatitude_1'].astype('float32')
  df['pointLongitude_2'] = df['pointLongitude_2'].astype('float32')
  df['pointLatitude_2'] = df['pointLatitude_2'].astype('float32')
  df['pointLongitude_3'] = df['pointLongitude_3'].astype('float32')
  df['pointLatitude_3'] = df['pointLatitude_3'].astype('float32')
  df['pointLongitude_4'] = df['pointLongitude_4'].astype('float32')
  df['pointLatitude_4'] = df['pointLatitude_4'].astype('float32')
  df['pointLongitude_5'] = df['pointLongitude_5'].astype('float32')
  df['pointLatitude_5'] = df['pointLatitude_5'].astype('float32')

  return df


def parseLine(line, outputType):

  returnValue = None
  if line.count(' - ') == 2:
    (dateStamp, logLevel, remainder) = line.split(' - ')
    if outputType == 'dateString':
      returnValue = datetime.strptime(dateStamp, 
        '%m/%d/%Y %H:%M:%S %p').isoformat() + 'Z'
    elif outputType == 'logLevel':
      returnValue = logLevel.split()
    if ':' in remainder:
      (key, value) = remainder.split(':', 1)
      if outputType == 'key':
        returnValue = key.strip()
      elif outputType == 'value':
        returnValue = value.strip()
  elif line.count(' - ') == 0 and line.count(': ') == 1:
    (key, value) = line.split(': ')
    if outputType == 'key':
      returnValue = key
    elif outputType == 'value':
      returnValue = value.strip()
  
  return returnValue


def getSentinelFiles(safe, meta):
  
  files = []
  doc = et.fromstring(meta['dataObjectSection'])
  dataObjects = doc.xpath('/dataObjectSection/dataObject')
  for ii in range(len(dataObjects)):
    file = {}
    href = doc.xpath('/dataObjectSection/dataObject[{0}]/byteStream/' \
      'fileLocation/@href'.format(ii+1))[0]
    fileParts = os.path.basename(href).split('-')
    file['name'] = os.path.join(safe, href[2:])
    file['ID'] = doc.xpath('/dataObjectSection/dataObject[{0}]/@ID' \
      .format(ii+1))[0]
    if file['ID'] == 'mapoverlay':
      file['type'] = 'preview'
    elif file['ID'] == 'productpreview':
      file['type'] = 'preview'
    elif file['ID'] == 'quicklook':
      file['type'] = 'preview'
    elif file['ID'].startswith('product'):
      file['type'] = 'annotation'
      file['swath'] = fileParts[1].upper()
      file['polarization'] = fileParts[3].upper()
      file['startTime'] = fileParts[4]
      file['stopTime'] = fileParts[5]
    elif file['ID'].startswith('calibration'):
      file['type'] = 'calibration'
      file['swath'] = fileParts[2].upper()
      file['polarization'] = fileParts[4].upper()
      file['startTime'] = fileParts[5]
      file['stopTime'] = fileParts[6]
    elif file['ID'].startswith('noise'):
      file['type'] = 'noise'
      file['swath'] = fileParts[2].upper()
      file['polarization'] = fileParts[4].upper()
      file['startTime'] = fileParts[5]
      file['stopTime'] = fileParts[6]
    elif file['ID'].startswith('rfi'):
      file['type'] = 'rfi'
      file['swath'] = fileParts[2].upper()
      file['polarization'] = fileParts[4].upper()
      file['startTime'] = fileParts[5]
      file['stopTime'] = fileParts[6]
    else:
      file['type'] = 'measurement'
      file['swath'] = fileParts[1].upper()
      file['polarization'] = fileParts[3].upper()
      file['startTime'] = fileParts[4]
      file['stopTime'] = fileParts[5]
    file['size'] = int(doc.xpath('/dataObjectSection/dataObject[{0}]/' \
      'byteStream/@size'.format(ii+1))[0])
    file['checksum'] = doc.xpath('/dataObjectSection/dataObject[{0}]/' \
      'byteStream/checksum'.format(ii+1))[0].text
    files.append(file)

  return files
