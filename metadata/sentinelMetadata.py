#!/usr/bin/python3

import zipfile
from datetime import datetime
import lxml.etree as et
from osgeo import gdal, osr


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
  