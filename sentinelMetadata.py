#!/usr/bin/python3

import os
import zipfile
from datetime import datetime
import lxml.etree as et
from osgeo import gdal, osr
import numpy as np
import glob


xsi = '{http://www.w3.org/2001/XMLSchema-instance}'
xfdu = '{urn:ccsds:schema:xfdu:1}'
gml = '{http://www.opengis.net/gml}'
safe = '{http://www.esa.int/safe/sentinel-1.0}'
s1 = '{http://www.esa.int/safe/sentinel-1.0/sentinel-1}'
s1sar = '{http://www.esa.int/safe/sentinel-1.0/sentinel-1/sar}'
s1sarl1 = '{http://www.esa.int/safe/sentinel-1.0/sentinel-1/sar/level-1}'
s1sarl2 = '{http://www.esa.int/safe/sentinel-1.0/sentinel-1/sar/level-2}'
gx = '{http://www.google.com/kml/ext/2.2}'
ns_xsi = {'xsi': 'http://www.w3.org/2001/XMLSchema-instance'}
ns_gml = {'gml': 'http://www.opengis.net/gml'}
ns_xfdu = {'xfdu': 'urn:ccsds:schema:xfdu:1'}
ns_safe = {'safe': 'http://www.esa.int/safe/sentinel-1.0'}
ns_s1 = {'s1': 'http://www.esa.int/safe/sentinel-1.0/sentinel-1'}
ns_s1sar = {'s1sar': 'http://www.esa.int/safe/sentinel-1.0/sentinel-1/sar'}
ns_s1sarl1 = \
  {'s1sarl1': 'http://www.esa.int/safe/sentinel-1.0/sentinel-1/sar/level-1'}
ns_s1sarl2 = \
  {'s1sarl2': 'http://www.esa.int/safe/sentinel-1.0/sentinel-1/sar/level-2'}
ns_gx = {'gx': 'http://www.google.com/kml/ext/2.2'}
ns = dict(
  list(ns_xsi.items()) +
  list(ns_gml.items()) +
  list(ns_xfdu.items()) +
  list(ns_safe.items()) +
  list(ns_s1.items()) +
  list(ns_s1sar.items()) +
  list(ns_s1sarl1.items()) +
  list(ns_s1sarl2.items()) +
  list(ns_gx.items())
)


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


def readManifestFile(manifestFile, zf=None):

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
  

def gammaRTClog2meta(logFile):

  ### Determine parent directory of log file
  parentDir = os.path.dirname(os.path.dirname(os.path.abspath(logFile)))

  ### Read information from log file
  fp = open(logFile, 'r')
  lines = fp.readlines()
  fp.close()

  meta = {}
  meta['ISO_RTC_metadataCreationTime'] = datetime.utcnow().isoformat() + 'Z'
  meta['ISO_RTC_productCreationTime'] = parseLine(lines[0], 'dateString')
  meta['ISO_RTC_inputGranule'] = \
    os.path.basename(parseLine(lines[1], 'value'))[:-5]
  meta['ISO_RTC_radiometry'] = parseLine(lines[3], 'value')
  meta['ISO_RTC_scale'] = parseLine(lines[4], 'value')
  outputBase = parseLine(lines[12], 'value')
  outputDir = os.path.join(parentDir, outputBase)
  meta['ISO_RTC_DEM_type'] = parseLine(lines[13], 'value').upper()
  for line in lines:
    value = parseLine(line, 'value')
    if value and value.startswith('geoid file'):
      meta['ISO_RTC_GammaVersion'] = \
        [i for i in value.split('/') if 'GAMMA_SOFTWARE' in i][0].split('-')[1]
    elif value and value.startswith('number of azimuth looks'):
      meta['ISO_RTC_azimuthLooks'] = int(value.split(':')[1])
    elif value and value.startswith('number of range looks'):
      meta['ISO_RTC_rangeLooks'] = int(value.split(':')[1])
  meta['ISO_RTC_XML_filename'] = outputBase + '.iso.xml'
    
  ### Extract more metadata from manifest file
  manifestFile = os.path.join(parentDir, 
    meta['ISO_RTC_inputGranule'] + '.SAFE', 'manifest.safe')
  manifest = readManifestFile(manifestFile)
  doc = et.fromstring(manifest['metadataSection'])
  meta['ISO_RTC_startTime'] = doc.xpath('/metadataSection/metadataObject' \
    '[@ID="acquisitionPeriod"]/metadataWrap/xmlData/' \
    'safe:acquisitionPeriod/safe:startTime', namespaces=ns)[0].text + 'Z'
  meta['ISO_RTC_stopTime'] = doc.xpath('/metadataSection/metadataObject' \
    '[@ID="acquisitionPeriod"]/metadataWrap/xmlData/' \
    'safe:acquisitionPeriod/safe:stopTime', namespaces=ns)[0].text + 'Z'
  
  ### Extract more metadata from annotation file
  annotationFile = glob.glob('{0}'.format(os.path.join(parentDir, 
    meta['ISO_RTC_inputGranule'] + '.SAFE', 'annotation', '*.xml')))[0]
  annotation = readAnnotationFile(annotationFile)
  doc = et.fromstring(annotation['adsHeader'])
  meta['ISO_RTC_mission'] = \
    doc.xpath('/adsHeader/missionId')[0].text.replace('S','SENTINEL-')
  meta['ISO_RTC_beamMode'] = doc.xpath('/adsHeader/mode')[0].text
  meta['ISO_RTC_absoluteOrbitNumber'] = \
    int(doc.xpath('/adsHeader/absoluteOrbitNumber')[0].text)
  doc = et.fromstring(annotation['generalAnnotation'])
  meta['ISO_RTC_orbitPassDirection'] = \
    doc.xpath('/generalAnnotation/productInformation/pass')[0].text.lower()
  speedOfLight = 299792458.0
  frequency = float(doc.xpath('/generalAnnotation/productInformation/' \
    'radarFrequency')[0].text)
  meta['ISO_RTC_wavelength'] = speedOfLight / frequency

  ### Extract Hyp3 version out of README file
  readmeFile = \
    os.path.join(parentDir, outputBase, outputBase + '.README.md.txt')
  with open(readmeFile, 'r') as fp:
    lines = fp.readlines()
  for line in lines:
    if parseLine(line, 'key') == 'Metadata version':
      meta['ISO_RTC_productVersion'] = parseLine(line, 'value')

  ### Setting additional metadata
  meta['ISO_RTC_lookDirection'] = 'RIGHT'
  if meta['ISO_RTC_mission'] == 'ALOS':
    meta['ISO_RTC_sensor'] = 'PALSAR'
  else:
    meta['ISO_RTC_sensor'] = 'SAR'

  ### Examine product directory for various files and calculate stats if exists
  polarization = meta['ISO_RTC_inputGranule'][14:16]
  if polarization == 'SH' or polarization == 'DH':
    mainPol = 'HH'
    polString = 'horizontal'
  elif polarization == 'SV' or polarization == 'DV':
    mainPol = 'VV'
    polString = 'vertical'

  ## Look at product file
  productFile = os.path.join(outputDir, outputBase + '_{0}.tif'.format(mainPol))
  if os.path.isfile(productFile):
    meta['ISO_RTC_terrainCorrectedImage'] = os.path.basename(productFile)
    raster = gdal.Open(productFile)
    meta['ISO_RTC_azimuthCount'] = raster.RasterXSize
    band = raster.GetRasterBand(1)
    stats = band.GetStatistics(False, True)
    proj = osr.SpatialReference()
    proj.ImportFromWkt(raster.GetProjectionRef())
    if proj.GetAttrValue('AUTHORITY', 0) == 'EPSG':
      epsg = int(proj.GetAttrValue('AUTHORITY', 1))
      geoTransform = raster.GetGeoTransform()
      meta['ISO_RTC_azimuthPixelSize'] = geoTransform[1]
      meta['ISO_RTC_rangePixelSize'] = -geoTransform[5]
      meta['ISO_RTC_azimuthCount'] = raster.RasterXSize
      meta['ISO_RTC_rangeCount'] = raster.RasterYSize
    else:
      epsg = 0
    meta['ISO_RTC_minValue'] = stats[0]
    meta['ISO_RTC_maxValue'] = stats[1]
    meta['ISO_RTC_meanValue'] = stats[2]
    meta['ISO_RTC_standardDeviation'] = stats[3]
    meta['ISO_RTC_transmittedPolarization'] = polString
    meta['ISO_RTC_receivedPolarization'] = polString
    meta['ISO_RTC_epsgCode'] = epsg

  ## Look at DEM file
  demFile = os.path.join(outputDir, outputBase + '_dem.tif')
  if os.path.isfile(demFile):
    meta['ISO_RTC_digitalElevationModel'] = os.path.basename(demFile)
    raster = gdal.Open(productFile)
    band = raster.GetRasterBand(1)
    stats = band.GetStatistics(False, True)
    meta['ISO_RTC_DEM_minValue'] = stats[0]
    meta['ISO_RTC_DEM_maxValue'] = stats[1]
    meta['ISO_RTC_DEM_meanValue'] = stats[2]
    meta['ISO_RTC_DEM_standardDeviation'] = stats[3]

  ## Look at the incidence angle map
  incidenceAngleFile = os.path.join(outputDir, outputBase + '_inc.tif')
  if os.path.isfile(incidenceAngleFile):
    meta['ISO_RTC_incidenceAngleMap'] = os.path.basename(incidenceAngleFile)
    raster = gdal.Open(productFile)
    band = raster.GetRasterBand(1)
    stats = band.GetStatistics(False, True)
    meta['ISO_RTC_incidenceMinValue'] = stats[0]
    meta['ISO_RTC_incidenceMaxValue'] = stats[1]
    meta['ISO_RTC_incidenceMeanValue'] = stats[2]
    meta['ISO_RTC_incidenceStandardDeviation'] = stats[3]

  ## Look at the layover shadow mask
  layoverShadowMaskFile = os.path.join(outputDir, outputBase + '_ls_map.tif')
  if os.path.isfile(layoverShadowMaskFile):
    meta['ISO_RTC_layoverShadowMask'] = os.path.basename(layoverShadowMaskFile)
    raster = gdal.Open(layoverShadowMaskFile)
    band = raster.GetRasterBand(1).ReadAsArray()
    pixelCount = raster.RasterXSize * raster.RasterYSize
    (histogram, _) = np.histogram(band, bins=np.arange(256))
    noLayoverShadow = 0.0
    trueLayover = 0.0
    layover = 0.0
    trueShadow = 0.0
    shadow = 0.0
    pixelCount -= float(histogram[0])
    for ii in range(len(histogram)):
      if ii & 1:
        noLayoverShadow += float(histogram[ii]) / pixelCount
      if ii & 2:
        trueLayover += float(histogram[ii]) / pixelCount
      if ii & 4:
        layover += float(histogram[ii]) / pixelCount
      if ii & 8:
        trueShadow += float(histogram[ii]) / pixelCount
      if ii & 16:
        shadow += float(histogram[ii]) / pixelCount
    meta['ISO_RTC_noLayoverShadowPercentage'] = noLayoverShadow * 100.0
    meta['ISO_RTC_trueLayoverPercentage'] = trueLayover * 100.0
    meta['ISO_RTC_layoverPercentage'] = layover * 100.0
    meta['ISO_RTC_trueShadowPercentage'] = trueShadow * 100.0
    meta['ISO_RTC_shadowPercentage'] = shadow * 100.0

  ## Look at the scattering area map
  scatteringAreaFile = os.path.join(outputDir, outputBase + '_area.tif')
  if os.path.isfile(scatteringAreaFile):
    meta['ISO_RTC_scatteringAreaMap'] = os.path.basename(scatteringAreaFile)
    raster = gdal.Open(productFile)
    band = raster.GetRasterBand(1)
    stats = band.GetStatistics(False, True)
    meta['ISO_RTC_scatteringMinValue'] = stats[0]
    meta['ISO_RTC_scatteringMaxValue'] = stats[1]
    meta['ISO_RTC_scatteringMeanValue'] = stats[2]
    meta['ISO_RTC_scatteringStandardDeviation'] = stats[3]

  ## Browse
  if polarization == 'DH' or polarization == 'DV':
    meta['ISO_RTC_browseImage'] = outputBase + '_rgb.png'
    meta['ISO_RTC_KMLoverlay'] = outputBase + '_rgb.kmz'
  else:
    meta['ISO_RTC_browseImage'] = outputBase + '.png'
    meta['ISO_RTC_KMLoverlay'] = outputBase + '.kmz'

  ## Extent
  (lon_min, lon_max, lat_min, lat_max) = get_latlon_extent(productFile)
  meta['ISO_RTC_westBoundLongitude'] = lon_min
  meta['ISO_RTC_eastBoundLongitude'] = lon_max
  meta['ISO_RTC_northBoundLatitude'] = lat_max
  meta['ISO_RTC_southBoundLatitude'] = lat_min

  return meta
