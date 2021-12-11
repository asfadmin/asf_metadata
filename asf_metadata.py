#!/usr/bin/python3

import os
import sys
from datetime import datetime, timedelta
import lxml.etree as et
from lxml import objectify
import json
import pandas as pd
import numpy as np
import configparser
import collections
import re
from osgeo import gdal, osr
from util import getValue, meta2dict, rreplace, uniqueList, \
  getParamsDataframe, params2dictList, upgradeDictionary2level, \
  mergeDictionaryParams, getLevelParamsList


def meta_project2epsg(metaDict):
  projection = metaDict['projection']['type']
  if projection == 'UNIVERSAL_TRANSVERSE_MERCATOR':
    hemisphere = metaDict['projection']['hem']
    zone = int(metaDict['projection']['param']['utm']['zone'])
    if hemisphere == 'N':
      epsg = 32600 + zone
    else:
      epsg = 32700 + zone
  return epsg


def update_template(template, listParams, listParamCounts):

  count = len(template)
  for listParam in listParams:
    indices = [i for i, x in enumerate(template) if x == listParam]
    indices.append(count)
    for ii in range(len(indices)-1):
      for kk in range(indices[ii], indices[ii+1]):
        if listParam in template[kk]:
          newElement = ('{0}[{1}]'.format(listParam, ii+1))
          template[kk] = rreplace(template[kk], listParam, newElement)

  return template


def iso_template2lists(excelFile, workSheet):

  ### ISO metadata structure
  meta = pd.read_excel(excelFile, sheet_name=workSheet)
  (rowCount, levelCount) = meta.shape
  metaType = list(meta['Type'])
  del meta['Type']
  metaValue = meta['Value']
  del meta['Value']
  if workSheet == 'ISO Metadata Structure':
    isoAttributes = list(meta['Attribute'])
    del meta['Attribute']
    isoAttributeValues = list(meta['AttributeValue'])
    del meta['AttributeValue']
  else:
    isoAttributes = []
    for ii in range(rowCount):
      isoAttributes.append('nan')
  (rowCount, levelCount) = meta.shape

  ### Read all entries in the template
  isoTemplate = []
  for ii in range(rowCount):
    param = ''
    for kk in range(levelCount):
      element = str(meta.at[ii, kk])
      if element != 'nan':
        param += element
        if kk < (levelCount-1) and str(meta.at[ii, kk+1]) != 'nan':
          param += '/'
    isoTemplate.append(param)

  ### Collect elements that contain lists
  duplicates = []
  for ii in range(levelCount):
    for kk in range(rowCount):
      occurences = [x.start() for x in re.finditer('/', isoTemplate[kk])]
      if len(occurences) == ii:
        if 'list' in metaType[kk]:
          duplicates.append(isoTemplate[kk])
    listParams = [item for item, \
      count in collections.Counter(duplicates).items() if count > 0]
    listParamCounts = {}
    for param in [element for index, element in enumerate(duplicates,1) 
      if element not in duplicates[index:]]:
        listParamCounts[param] = duplicates.count(param)
    if len(listParams) > 0:
      isoTemplate = update_template(isoTemplate, listParams, listParamCounts)

  ### Read parameters out of template
  isoParams = []
  isoValues = []
  for ii in range(rowCount):
    if 'attribute' in metaType[ii]:
      attrib = isoTemplate[ii] + '/' + str(isoAttributes[ii])
      isoParams.append(attrib)
      isoValues.append(isoAttributeValues[ii])
    if 'value' in metaType[ii]:
      isoParams.append(isoTemplate[ii])
      isoValues.append(metaValue[ii])

  return (isoTemplate, isoParams, isoValues)


def merge_iso_dem_template(isoTemplate, demTemplate):

  ### Analyze ISO template
  maxIdentificationInfo = 0
  maxIdentificationIndex = 0
  maxContentInfo = 0
  maxContentIndex = 0
  for ii in range(len(isoTemplate)):
    testLevel = isoTemplate[ii][-2:-1]
    if isoTemplate[ii][:-3].endswith('identificationInfo'):
      if maxIdentificationInfo < int(testLevel):
        maxIdentificationInfo = int(testLevel)
    if 'identificationInfo' in isoTemplate[ii] and maxIdentificationIndex < ii:
      maxIdentificationIndex = ii
    if isoTemplate[ii][:-3].endswith('contentInfo'):
      if maxContentInfo < int(testLevel):
        maxContentInfo = int(testLevel)
    if 'contentInfo' in isoTemplate[ii] and maxContentIndex < ii:
      maxContentIndex = ii

  ### Analyze DEM template
  maxIndex = 0
  for ii in range(len(demTemplate)):
    if 'identificationInfo' in demTemplate[ii] and maxIndex < ii:
      maxIndex = ii

  ### Combine ISO and DEM templates
  mergedTemplate = []
  for ii in range(maxIdentificationIndex+1):
    mergedTemplate.append(isoTemplate[ii])
  originalString = 'identificationInfo[1]'
  replaceString = ('identificationInfo[{0}]'.format(maxIdentificationInfo+1))
  for ii in range(maxIndex+1):
    if 'identificationInfo' in demTemplate[ii]:
      demInfo = demTemplate[ii].replace(originalString, replaceString)
      mergedTemplate.append(demInfo)
  for ii in range(maxIdentificationIndex+1, maxContentIndex+1):
    mergedTemplate.append(isoTemplate[ii])
  originalString = 'contentInfo[1]'
  replaceString = ('contentInfo[{0}]'.format(maxIdentificationInfo+1))
  for ii in range(maxIndex+1, len(demTemplate)):
    demInfo = demTemplate[ii].replace(originalString, replaceString)
    mergedTemplate.append(demInfo)
  for ii in range(maxContentIndex+1, len(isoTemplate)):
    mergedTemplate.append(isoTemplate[ii])

  return mergedTemplate


def add_dem_lists(isoTemplate, isoParams, isoValues, 
  demTemplate, demParams, demValues):

  ### Analyze ISO template
  maxIdentificationInfo = 0
  maxIdentificationIndex = 0
  maxContentInfo = 0
  maxContentIndex = 0
  for ii in range(len(isoTemplate)):
    testLevel = isoTemplate[ii][-2:-1]
    if isoTemplate[ii][:-3].endswith('identificationInfo'):
      if maxIdentificationInfo < int(testLevel):
        maxIdentificationInfo = int(testLevel)
    if 'identificationInfo' in isoTemplate[ii] and maxIdentificationIndex < ii:
      maxIdentificationIndex = ii
    if isoTemplate[ii][:-3].endswith('contentInfo'):
      if maxContentInfo < int(testLevel):
        maxContentInfo = int(testLevel)
    if 'contentInfo' in isoTemplate[ii] and maxContentIndex < ii:
      maxContentIndex = ii

  ### Analyze DEM template
  maxIndex = 0
  for ii in range(len(demTemplate)):
    if 'identificationInfo' in demTemplate[ii] and maxIndex < ii:
      maxIndex = ii

  ### Combine ISO and DEM templates
  mergedTemplate = []
  for ii in range(maxIdentificationIndex+1):
    mergedTemplate.append(isoTemplate[ii])
  originalString = 'identificationInfo[1]'
  replaceString = ('identificationInfo[{0}]'.format(maxIdentificationInfo+1))
  for ii in range(maxIndex+1):
    if 'identificationInfo' in demTemplate[ii]:
      demInfo = demTemplate[ii].replace(originalString, replaceString)
      mergedTemplate.append(demInfo)
  for ii in range(maxIdentificationIndex+1, maxContentIndex+1):
    mergedTemplate.append(isoTemplate[ii])
  originalString = 'contentInfo[1]'
  replaceString = ('contentInfo[{0}]'.format(maxIdentificationInfo+1))
  for ii in range(maxIndex+1, len(demTemplate)):
    demInfo = demTemplate[ii].replace(originalString, replaceString)
    mergedTemplate.append(demInfo)
  for ii in range(maxContentIndex+1, len(isoTemplate)):
    mergedTemplate.append(isoTemplate[ii])

  ### Analyze ISO params and values
  maxIdentificationIndex = 0
  maxContentIndex = 0
  for ii in range(len(isoParams)):
    if 'identificationInfo' in isoParams[ii] and maxIdentificationIndex < ii:
      maxIdentificationIndex = ii
    if 'contentInfo' in isoParams[ii] and maxContentIndex < ii:
      maxContentIndex = ii

  ### Analyze DEM params
  maxIndex = 0
  for ii in range(len(demParams)):
    if 'identificationInfo' in demParams[ii] and maxIndex < ii:
      maxIndex = ii

  ### Combine ISO and DEM params and values
  mergedParams = []
  mergedValues = []
  for ii in range(maxIdentificationIndex+1):
    mergedParams.append(isoParams[ii])
    mergedValues.append(isoValues[ii])
  originalString = 'identificationInfo[1]'
  replaceString = ('identificationInfo[{0}]'.format(maxIdentificationInfo+1))
  for ii in range(maxIndex+1):
    demInfo = demParams[ii].replace(originalString, replaceString)
    mergedParams.append(demInfo)
    mergedValues.append(demValues[ii])
  for ii in range(maxIdentificationIndex+1, maxContentIndex+1):
    mergedParams.append(isoParams[ii])
    mergedValues.append(isoValues[ii])
  originalString = 'contentInfo[1]'
  replaceString = ('contentInfo[{0}]'.format(maxIdentificationInfo+1))
  for ii in range(maxIndex+1, len(demParams)):
    demInfo = demParams[ii].replace(originalString, replaceString)
    mergedParams.append(demInfo)
    mergedValues.append(demValues[ii])
  for ii in range(maxContentIndex+1, len(isoParams)):
    mergedParams.append(isoParams[ii])
    mergedValues.append(isoValues[ii])

  return (mergedTemplate, mergedParams, mergedValues)


def iso_attributes2lists(excelFile):

  ### ISO metadata structure
  meta = pd.read_excel(excelFile, sheet_name='ISO Metadata Structure')
  isoAttributes = list(meta['Attribute'])
  isoAttributeValues = list(meta['AttributeValue'])
  (isoTemplate, _, _) = iso_template2lists(excelFile, 'ISO Metadata Structure')

  return (isoAttributes, isoAttributeValues, isoTemplate)


def iso_xml_structure(excelFile, isoTemplate, isoParams, isoValues, nsFlag):

  ### Get namespace look up table
  nsLUT = pd.read_excel(excelFile, sheet_name='Namespaces LUT')
  nsPairs = dict(zip(nsLUT['Parameter'], nsLUT['Namespace']))

  ### Build namespace
  nsLUT = pd.read_excel(excelFile, sheet_name='Namespaces URL')
  nsLUT = nsLUT.set_index('Namespace')
  ns_xsi = {'xsi': nsLUT.at['xsi','URL']}
  ns_gmd = {'gmd': nsLUT.at['gmd','URL']}
  ns_gco = {'gco': nsLUT.at['gco','URL']}
  ns_xs = {'xs': nsLUT.at['xs','URL']}
  ns_eos = {'eos': nsLUT.at['eos','URL']}
  ns_echo = {'echo': nsLUT.at['echo','URL']}
  ns_xlink = {'xlink': nsLUT.at['xlink','URL']}
  ns_gml = {'gml': nsLUT.at['gml','URL']}
  ns_gmi = {'gmi': nsLUT.at['gmi','URL']}
  ns_gmx = {'gmx': nsLUT.at['gmx','URL']}
  ns = dict(
    list(ns_xsi.items()) +
    list(ns_gmd.items()) +
    list(ns_gco.items()) +
    list(ns_xs.items()) +
    list(ns_eos.items()) +
    list(ns_echo.items()) +
    list(ns_xlink.items()) +
    list(ns_gml.items()) +
    list(ns_gmi.items()) +
    list(ns_gmx.items())
  )

  isoElements = []
  isoElementCount = []
  rowCount = len(isoTemplate)
  levelCount = 0
  for ii in range(rowCount):
    metaElements = isoTemplate[ii].split('/')
    if len(metaElements) > levelCount:
      levelCount = len(metaElements)
    isoElements.append(metaElements)
    isoElementCount.append(len(metaElements))
  levelElements = []
  for kk in range(levelCount):
    meta = []
    for ii in range(rowCount):
      if kk < len(isoElements[ii]):
        meta.append(isoElements[ii][kk])
    levelElements.append(uniqueList(meta))

  ### Build XMl element tree from nested dictionary
  nameSpace = ('{%s}' % nsLUT.at[nsPairs['DS_Series'],'URL'])
  if nsFlag == True:
    root = et.Element('{0}{1}'.format(nameSpace, 'DS_Series'), nsmap=ns)
  else:
    root = et.Element('DS_Series')

  ## Build structure
  ref = {}
  ref['DS_Series'] = root
  for ii in range(1, levelCount):
    indices = [kk for kk, x in enumerate(isoElementCount) if x == ii+1]
    for kk in range(len(indices)):
      element = isoTemplate[indices[kk]]
      parent = ref[element.rsplit('/',1)[0]]
      param = element.rsplit('/',1)[1]
      if 'EOS_AdditionalAttributes' in element:
        nameSpace = ('{%s}' % nsLUT.at['eos','URL'])
      elif ('acquisitionInformation' in element) and \
        ((param == 'identifier') or (param == 'description') or \
        (param == 'type')):
        nameSpace = ('{%s}' % nsLUT.at['gmi','URL'])
      else:
        nameSpace = ('{%s}' % nsLUT.at[nsPairs[param.split('[')[0]],'URL'])
      if nsFlag == True:
        child = et.SubElement(parent, '{0}{1}'.format(nameSpace,
          param.split('[')[0], nsmap=ns))
      else:
        child = et.SubElement(parent, param.split('[')[0])
      ref[element] = child

  ## Add values
  for ii in range(len(isoParams)):
    element = isoParams[ii]
    param = element.rsplit('/',1)[1]
    if (param == 'value') and (('result/value' not in element) or \
      ('EOS_AdditionalAttribute/value' not in element)):
      continue
    elif param == 'id':
      reference = ref[element.rsplit('/',1)[0]]
      nameSpace = ('{%s}' % nsLUT.at[nsPairs[param.split('[')[0]],'URL'])
      if nsFlag == True:
        nsParam = '{0}{1}'.format(nameSpace, param.split('[')[0], nsmap=ns)
      else:
        nsParam = param.split('[')[0]
      reference.attrib[nsParam] = str(isoValues[ii])
    elif (param == 'href') or (param == 'nilReason'):
      reference = ref[element.rsplit('/',1)[0]]
      nameSpace = ('{%s}' % nsLUT.at[nsPairs[param.split('[')[0]],'URL'])
      if nsFlag == True:
        nsParam = '{0}{1}'.format(nameSpace, param.split('[')[0], nsmap=ns)
      else:
        nsParam = param.split('[')[0]
      reference.attrib[nsParam] = str(isoValues[ii])
    elif param == 'type':
      reference = ref[element.rsplit('/',1)[0]]
      nameSpace = ('{%s}' % nsLUT.at['xsi','URL'])
      if nsFlag == True:
        nsParam = '{0}{1}'.format(nameSpace, param.split('[')[0], nsmap=ns)
      else:
        nsParam = param.split('[')[0]
      reference.attrib[nsParam] = str(isoValues[ii])
    elif param == 'uom':
      reference = ref[element.rsplit('/',1)[0]]
      reference.attrib[param] = str(isoValues[ii])
      reference.text = str(isoValues[ii-1])
    elif param.split('[')[0] in nsPairs:
      nameSpace = ('{%s}' % nsLUT.at[nsPairs[param.split('[')[0]],'URL'])
      reference = ref[element]
      reference.text = str(isoValues[ii])
      if 'Code' in param:
        if 'EOS' in element:
          codeURL = 'http://earthdata.nasa.gov/metadata/resources/' \
            'Codelist.xml#' + param
        else:
          codeURL = 'https://cdn.earthdata.nasa.gov/iso/resources/Codelist/' \
            'gmxCodelists.xml#' + param
        reference.attrib['codeList'] = codeURL
        reference.attrib['codeListValue'] = str(isoValues[ii])
      if param == 'RecordType':
        nameSpace = ('{%s}' % nsLUT.at['xlink','URL'])
        if nsFlag == True:
          nsParam = '{0}{1}'.format(nameSpace, 'href', nsmap=ns)
        else:
          nsParam = param.split('[')[0]
        reference.attrib[nsParam] = "http://earthdata.nasa.gov/schemas/eos/" \
          "eos.xsd#xpointer(//element[@name='EOS_AdditionalAttributes'])"

  return root


def add_dem_attributes(isoAttributes, demAttributes, isoTemplate, demTemplate):

  ### Analyze ISO template
  maxIdentificationInfo = 0
  maxIdentificationIndex = 0
  maxContentInfo = 0
  maxContentIndex = 0
  for ii in range(len(isoTemplate)):
    testLevel = isoTemplate[ii][-2:-1]
    if isoTemplate[ii][:-3].endswith('identificationInfo'):
      if maxIdentificationInfo < int(testLevel):
        maxIdentificationInfo = int(testLevel)
    if 'identificationInfo' in isoTemplate[ii] and maxIdentificationIndex < ii:
      maxIdentificationIndex = ii
    if isoTemplate[ii][:-3].endswith('contentInfo'):
      if maxContentInfo < int(testLevel):
        maxContentInfo = int(testLevel)
    if 'contentInfo' in isoTemplate[ii] and maxContentIndex < ii:
      maxContentIndex = ii

  ### Analyze DEM template
  maxIndex = 0
  for ii in range(len(demTemplate)):
    if 'identificationInfo' in demTemplate[ii] and maxIndex < ii:
      maxIndex = ii

  ### Combine ISO and DEM attributes
  mergedAttributes = []
  for ii in range(maxIdentificationIndex+1):
    mergedAttributes.append(isoAttributes[ii])
  for ii in range(maxIndex+1):
    mergedAttributes.append(demAttributes[ii])
  for ii in range(maxIdentificationIndex+1, maxContentIndex+1):
    mergedAttributes.append(isoAttributes[ii])
  for ii in range(maxIndex+1, len(demAttributes)):
    mergedAttributes.append(demAttributes[ii])
  for ii in range(maxContentIndex+1, len(isoAttributes)):
    mergedAttributes.append(isoAttributes[ii])

  return mergedAttributes


def iso_dictionary_structure(excelFile, demFile, metaTemplate, metaParams, 
  metaValues):

  ### Read attributes and their values from Excel spreadsheet
  (isoAttributes, _, isoTemplate) = iso_attributes2lists(excelFile)
  if demFile:
    (demAttributes, _, demTemplate) = iso_attributes2lists(demFile)
    isoAttributes = add_dem_attributes(isoAttributes, demAttributes[5:], 
      isoTemplate, demTemplate[5:])

  ### Build dataframe with parameters
  dfParams = getParamsDataframe(metaParams)
  (_, oldLevel) = dfParams.shape
  oldLevel -= 1
  dictList = \
    params2dictList(dfParams, metaParams, metaValues, oldLevel)

  ### Build dictionary structure - one level at a time
  for newLevel in reversed(range(oldLevel)):

    newDictList = \
      params2dictList(dfParams, metaParams, metaValues, newLevel)
    dictList += newDictList
    dictList = upgradeDictionary2level(dictList, oldLevel, newLevel)
    levelParams = \
      getLevelParamsList(metaTemplate, isoAttributes, newLevel)
    dictList = mergeDictionaryParams(dictList, levelParams, oldLevel, newLevel)
    oldLevel = newLevel

  isoDict = {}
  isoDict[dictList[0]['key']] = dictList[0]['dict']

  return isoDict


def umm_dictionary_structure(metaTemplate, metaParams, metaValues):

  ### Build dataframe with parameters
  dfParams = getParamsDataframe(metaParams)
  (_, oldLevel) = dfParams.shape
  oldLevel -= 1
  dictList = \
    params2dictList(dfParams, metaParams, metaValues, oldLevel)

  ### Build dictionary structure - one level at a time
  for newLevel in reversed(range(oldLevel)):

    newDictList = \
      params2dictList(dfParams, metaParams, metaValues, newLevel)
    dictList += newDictList
    dictList = upgradeDictionary2level(dictList, oldLevel, newLevel)
    levelParams = getLevelParamsList(metaTemplate, None, newLevel)
    dictList = mergeDictionaryParams(dictList, levelParams, oldLevel, newLevel)
    oldLevel = newLevel

  isoDict = {}
  isoDict[dictList[0]['key']] = dictList[0]['dict']

  return isoDict


def meta_json_file(metaStructure, jsonFile):

  ### Write JSON file
  with open(jsonFile, 'w') as outF:
    outF.write(json.dumps(metaStructure, indent=2))
  outF.close()


def meta_xml_file(metaStructure, xmlFile):

  ### Write XML file
  with open(xmlFile, 'wb') as outF:
    outF.write(et.tostring(metaStructure, xml_declaration=True, encoding='utf-8',
      pretty_print=True))
  outF.close()


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


def gamma_rtc_metadata(config):
	
  timestamp = datetime.utcnow().isoformat() + 'Z'

  ## Get relevant file names from processing list file
  fileName = {}
  fileName['granule'] = config['GAMMA RTC']['granule']
  fileName['metadata'] = config['GAMMA RTC']['metadata']
  fileName['oversampled dem file'] = \
    config['GAMMA RTC']['oversampled dem file']
  fileName['oversampled dem metadata'] = \
    config['GAMMA RTC']['oversampled dem metadata']
  fileName['original dem file'] = config['GAMMA RTC']['original dem file']
  fileName['layover shadow mask'] = config['GAMMA RTC']['layover shadow mask']
  fileName['layover shadow stats'] = \
    config['GAMMA RTC']['layover shadow stats']
  fileName['layover shadow stats'] = \
    config['GAMMA RTC']['layover shadow stats']
  fileName['incidence angle file'] = \
    config['GAMMA RTC']['incidence angle file']
  fileName['incidence angle metadata'] = \
    config['GAMMA RTC']['incidence angle metadata']
  if config.has_option('GAMMA RTC', 'input HH file') == True:
    fileName['input HH file'] = config['GAMMA RTC']['input HH file']
  if config.has_option('GAMMA RTC', 'input HV file') == True:
    fileName['input HV file'] = config['GAMMA RTC']['input HV file']
  if config.has_option('GAMMA RTC', 'input VH file') == True:
    fileName['input VH file'] = config['GAMMA RTC']['input VH file']
  if config.has_option('GAMMA RTC', 'input VV file') == True:
    fileName['input VV file'] = config['GAMMA RTC']['input VV file']
  if config.has_option('GAMMA RTC', 'terrain corrected HH metadata') == True:
    fileName['terrain corrected HH metadata'] = \
      config['GAMMA RTC']['terrain corrected HH metadata']
  if config.has_option('GAMMA RTC', 'terrain corrected HV metadata') == True:
    fileName['terrain corrected HV metadata'] = \
      config['GAMMA RTC']['terrain corrected VH metadata']
  if config.has_option('GAMMA RTC', 'terrain corrected VH metadata') == True:
    fileName['terrain corrected VH metadata'] = \
      config['GAMMA RTC']['terrain corrected VH metadata']
  if config.has_option('GAMMA RTC', 'terrain corrected VV metadata') == True:
    fileName['terrain corrected VV metadata'] = \
      config['GAMMA RTC']['terrain corrected VV metadata']
  if config.has_option('GAMMA RTC', 'terrain corrected HH file') == True:
    fileName['terrain corrected HH file'] = \
      config['GAMMA RTC']['terrain corrected HH file']
  if config.has_option('GAMMA RTC', 'terrain corrected HV file') == True:
    fileName['terrain corrected HV file'] = \
      config['GAMMA RTC']['terrain corrected HV file']
  if config.has_option('GAMMA RTC', 'terrain corrected VH file') == True:
    fileName['terrain corrected VH file'] = \
      config['GAMMA RTC']['terrain corrected VH file']
  if config.has_option('GAMMA RTC', 'terrain corrected VV file') == True:
    fileName['terrain corrected VV file'] = \
      config['GAMMA RTC']['terrain corrected VV file']
  fileName['initial processing log'] = \
    config['GAMMA RTC']['initial processing log']
  fileName['terrain correction log'] = \
    config['GAMMA RTC']['terrain correction log']
  fileName['main log'] = config['GAMMA RTC']['main log']
  fileName['mk_geo_radcal_0 log'] = config['GAMMA RTC']['mk_geo_radcal_0 log']
  fileName['mk_geo_radcal_1 log'] = config['GAMMA RTC']['mk_geo_radcal_1 log']
  fileName['mk_geo_radcal_2 log'] = config['GAMMA RTC']['mk_geo_radcal_2 log']
  fileName['mk_geo_radcal_3 log'] = config['GAMMA RTC']['mk_geo_radcal_3 log']
  fileName['coreg_check log'] = config['GAMMA RTC']['coreg_check log']
  fileName['mli.par file'] = config['GAMMA RTC']['mli.par file']
  fileName['browse image'] = config['GAMMA RTC']['browse image']
  fileName['kml overlay'] = config['GAMMA RTC']['kml overlay']
  fileName['gamma version'] = config['GAMMA RTC']['gamma version']
  fileName['hyp3_rtc version'] = config['GAMMA RTC']['hyp3_rtc version']
  demSource = config['GAMMA RTC']['dem source']

  ## Check file existence
  if os.path.isfile(fileName['metadata']) == False:
    print('Metadata file ({0}) is missing!'.format(fileName['metadata']))
    sys.exit(1)
  if os.path.isfile(fileName['oversampled dem file']) == False:
    print('Oversampled digital elevation model ({0}) is missing!' \
      .format(fileName['oversampled dem file']))
  if os.path.isfile(fileName['original dem file']) == False:
    print('Original digital elevation model ({0}) is missing!' \
      .format(fileName['original dem file']))
  if ('terrain corrected HH metadata' in fileName) and \
    (os.path.isfile(fileName['terrain corrected HH metadata']) == False):
    print('Metadata for terrain corrected file ({0}) is missing!' \
      .format(fileName['terrain corrected HH metadata']))
    sys.exit(1)
  if ('terrain corrected HH file' in fileName) and \
    (os.path.isfile(fileName['terrain corrected HH file']) == False):
    print('Terrain corrected file ({0}) is missing!' \
      .format(fileName['terrain corrected HH file']))
    sys.exit(1)
  if ('terrain corrected HV metadata' in fileName) and \
    (os.path.isfile(fileName['terrain corrected HV metadata']) == False):
    print('Metadata for terrain corrected file ({0}) is missing!' \
      .format(fileName['terrain corrected HV metadata']))
    sys.exit(1)
  if ('terrain corrected HV file' in fileName) and \
    (os.path.isfile(fileName['terrain corrected HV file']) == False):
    print('Terrain corrected file ({0}) is missing!' \
      .format(fileName['terrain corrected HV file']))
    sys.exit(1)
  if ('terrain corrected VH metadata' in fileName) and \
    (os.path.isfile(fileName['terrain corrected VH metadata']) == False):
    print('Metadata for terrain corrected file ({0}) is missing!' \
      .format(fileName['terrain corrected VH metadata']))
    sys.exit(1)
  if ('terrain corrected VH file' in fileName) and \
    (os.path.isfile(fileName['terrain corrected VH file']) == False):
    print('Terrain corrected file ({0}) is missing!' \
      .format(fileName['terrain corrected VH file']))
    sys.exit(1)
  if ('terrain corrected VV metadata' in fileName) and \
    (os.path.isfile(fileName['terrain corrected VV metadata']) == False):
    print('Metadata for terrain corrected file ({0}) is missing!' \
      .format(fileName['terrain corrected VV metadata']))
    sys.exit(1)
  if ('terrain corrected VV file' in fileName) and \
    (os.path.isfile(fileName['terrain corrected VV file']) == False):
    print('Terrain corrected file ({0}) is missing!' \
      .format(fileName['terrain corrected VV file']))
    sys.exit(1)
  if os.path.isfile(fileName['browse image']) == False:
    print('Browse image file ({0}) is missing!' \
      .format(fileName['browse image']))
  if os.path.isfile(fileName['kml overlay']) == False:
    print('KML overlay file ({0}) is missing!' \
      .format(fileName['kml overlay']))

  ## Determine the values for the ISO metadata
  meta = {}
  meta['ISO_RTC_outputFile'] = fileName['granule'] + '.tif'
  meta['ISO_RTC_XML_filename'] = fileName['granule'] + '.iso.xml'
  meta['ISO_RTC_metadataCreationTime'] = timestamp
  meta['ISO_RTC_digitalElevationModel'] = \
    os.path.basename(fileName['original dem file'])
  meta['ISO_RTC_layoverShadowMask'] = \
    os.path.basename(fileName['layover shadow mask'])
  meta['ISO_RTC_incidenceAngleMap'] = \
    os.path.basename(fileName['incidence angle file'])

  if os.path.isfile(fileName['mli.par file']) == True:
    print('Extracting information from {0} ...' \
      .format(fileName['mli.par file']))
    lines = [line.rstrip() for line in open(fileName['mli.par file'])]
    for line in lines:
      if ':' in line:
        param = line.split(':')[0]
        value = line.split(':')[1].lstrip()
        if param == 'sensor':
          if 'ERS1' in value:
            platform = 'ERS-1'
          elif 'ERS2' in value:
            platform = 'ERS-2'
          elif 'PALSAR' in value:
            platform = 'ALOS'
          elif 'S1A' in value:
            platform = 'Sentinel-1A'
          elif 'S1B' in value:
            platform = 'Sentinel-1B'
          else:
            print('Could not identify sensor!')
            sys.exit(1)
        if param == 'date':
          dateTimestamp = datetime.strptime(value, '%Y %m %d %H %M %S.%f')
          dateStamp = datetime.combine(dateTimestamp.date(),
            datetime.min.time())
        if param == 'start_time':
          startTime = dateStamp + timedelta(seconds=float(value.split()[0]))
        if param == 'end_time':
          endTime = dateStamp + timedelta(seconds=float(value.split()[0]))
        if param == 'radar_frequency':
          speedOfLight = 299792458.0
          frequency = float(value.split()[0])
          wavelength = speedOfLight / frequency
        if param == 'range_looks':
          rangeLooks = int(value)
        if param == 'azimuth_looks':
          azimuthLooks = int(value)
        if param == 'range_pixel_spacing':
          slantSpacing = float(value.split()[0])
        if param == 'azimuth_pixel_spacing':
          azimuthSpacing = float(value.split()[0])
  else:
    print('Could not find MLI parameter file ({0})!' \
      .format(fileName['mli.par file']))
    sys.exit(1)

  if os.path.isfile(fileName['mk_geo_radcal_2 log']) == True:
    print('Extracting information from {0} ...' \
      .format(fileName['mk_geo_radcal_2 log']))
    lines = [line.rstrip() for line in open(fileName['mk_geo_radcal_2 log'])]
    for line in lines:
      if ':' in line:
        param = line.split(':')[0]
        value = line.split(':',1)[1].lstrip()
        if param == 'final range offset poly. coeff.':
          rangeOffset = float(value.split()[0])
        if param == 'final azimuth offset poly. coeff.':
          azimuthOffset = float(value.split()[0])
        if param == 'final solution':
          patches = value.split()
  else:
    print('Could not find offset information ({0}) !' \
      .format(fileName['mk_geo_radcal_2 log']))
    sys.exit(1)

  coregistrationPassed = -1
  if os.path.isfile(fileName['coreg_check log']) == True:
    lines = [line.rstrip() for line in open(fileName['coreg_check log'])]
    for line in lines:
      if 'passed coregistration' in line:
        coregistrationPassed = True
      if 'failed coregistration' in line:
        coregistrationPassed = False
    if coregistrationPassed == -1:
      print('Could not determine if this granule passed coregistratioon!')
      sys.exit(1)
  else:
    print('Could not find coregistration check log file ({0}) !' \
      .format(fileName['coreg_check log']))

  # Metadata - Input image
  if os.path.isfile(fileName['metadata']) == True:
    metaDict = meta2dict(fileName['metadata'])
    if platform == 'ALOS':
      if 'input HH file' in fileName:
        originalFile = os.path.basename(fileName['input HH file'])
      elif 'input VV file' in fileName:
        originalFile = os.path.basename(fileName['input VV file'])
      frame = int(metaDict['general']['frame'])
    elif 'Sentinel' in platform:
      originalFile = metaDict['general']['name'].replace('.iso.xml','')
      orbit = int(metaDict['general']['name'][49:55])
      frame = -1
    else:
      print('Unsupported platform: {0}'.format(platform))
      sys.exit(1)

    meta['ISO_RTC_inputGranule'] = originalFile
    meta['ISO_RTC_mission'] = platform
    if platform == 'ALOS':
      meta['ISO_RTC_sensor'] = 'PALSAR'
    else:
      meta['ISO_RTC_sensor'] = 'SAR'
    meta['ISO_RTC_lookDirection'] = 'RIGHT'
    meta['ISO_RTC_wavelength'] = wavelength
    meta['ISO_RTC_beamMode'] = metaDict['general']['mode']
    if metaDict['general']['orbit_direction'] == 'A':
      meta['ISO_RTC_orbitDirection'] = 'ASCENDING'
    elif metaDict['general']['orbit_direction'] == 'D':
      meta['ISO_RTC_orbitDirection'] = 'DESCENDING'
    meta['ISO_RTC_absoluteOrbitNumber'] = orbit
    if frame > 0:
      meta['ISO_RTC_frame'] = frame
    if metaDict['general']['orbit_direction'] == 'A':
      meta['ISO_RTC_orbitPassDirection'] = 'ascending'
    elif metaDict['general']['orbit_direction'] == 'D':
      meta['ISO_RTC_orbitPassDirection'] = 'descending'
    meta['ISO_RTC_startTime'] = startTime.isoformat() + 'Z'
    meta['ISO_RTC_stopTime'] = endTime.isoformat() + 'Z'
    meta['ISO_RTC_rangeLooks'] = rangeLooks
    meta['ISO_RTC_azimuthLooks'] = azimuthLooks
    meta['ISO_RTC_rangePixelSize'] = slantSpacing
    meta['ISO_RTC_azimuthPixelSize'] = azimuthSpacing
    metaDict = None
  else:
    print('Metadata file ({0}) does not exist !'.format(fileName['metadata']))
    sys.exit(1)

  # Metadata - terrain correction
  if 'terrain corrected HH metadata' in fileName:
    rtcFilename = os.path.basename(fileName['terrain corrected HH file'])
    metaDict = meta2dict(fileName['terrain corrected HH metadata'])
  elif 'terrain corrected VV metadata' in fileName:
    rtcFilename = os.path.basename(fileName['terrain corrected VV file'])
    metaDict = meta2dict(fileName['terrain corrected VV metadata'])
  else:
    print('Could not find metadata for terrain corrected product !')
    sys.exit(1)

  meta['ISO_RTC_terrainCorrectedImage'] = rtcFilename
  meta['ISO_RTC_azimuthCount'] = metaDict['general']['sample_count']
  meta['ISO_RTC_azimuthPixelSize'] = metaDict['general']['y_pixel_size']
  meta['ISO_RTC_rangeCount'] = metaDict['general']['line_count']
  meta['ISO_RTC_rangePixelSize'] = metaDict['general']['x_pixel_size']
  meta['ISO_RTC_epsgCode'] = meta_project2epsg(metaDict)
  metaDict = None

  # Metadata - terrain correction
  meta['ISO_RTC_GammaVersion'] = fileName['gamma version']
  meta['ISO_RTC_productVersion'] = fileName['hyp3_rtc version']
  meta['ISO_RTC_coregistrationPatches'] = int(patches[6])
  meta['ISO_RTC_successfulPatches'] = int(patches[0])
  meta['ISO_RTC_coregistrationSuccessFlag'] = coregistrationPassed
  meta['ISO_RTC_coregistrationRangeOffset'] = rangeOffset
  meta['ISO_RTC_coregistrationAzimuthOffset'] = azimuthOffset

  # Browse
  meta['ISO_RTC_browseImage'] = os.path.basename(fileName['browse image'])
  meta['ISO_RTC_KMLoverlay'] = os.path.basename(fileName['kml overlay'])

  # Extent
  if 'terrain corrected HH file' in fileName:
    (lon_min, lon_max, lat_min, lat_max) = \
      get_latlon_extent(fileName['terrain corrected HH file'])
  else:
    (lon_min, lon_max, lat_min, lat_max) = \
      get_latlon_extent(fileName['terrain corrected VV file'])
  meta['ISO_RTC_westBoundLongitude'] = lon_min
  meta['ISO_RTC_eastBoundLongitude'] = lon_max
  meta['ISO_RTC_northBoundLatitude'] = lat_max
  meta['ISO_RTC_southBoundLatitude'] = lat_min

  # Statistics - only for the first full-pol RTC image
  if 'terrain corrected HH file' in fileName:
    raster = gdal.Open(fileName['terrain corrected HH file'])
    band = raster.GetRasterBand(1)
    stats = band.GetStatistics(False, True)
    meta['ISO_RTC_minValue'] = stats[0]
    meta['ISO_RTC_maxValue'] = stats[1]
    meta['ISO_RTC_meanValue'] = stats[2]
    meta['ISO_RTC_standardDeviation'] = stats[3]
    meta['ISO_RTC_transmittedPolarization'] = 'horizontal'
    meta['ISO_RTC_receivedPolarization'] = 'horizontal'
  elif 'terrain corrected VV file' in fileName:
    raster = gdal.Open(fileName['terrain corrected VV file'])
    band = raster.GetRasterBand(1)
    stats = band.GetStatistics(False, True)
    meta['ISO_RTC_minValue'] = stats[0]
    meta['ISO_RTC_maxValue'] = stats[1]
    meta['ISO_RTC_meanValue'] = stats[2]
    meta['ISO_RTC_standardDeviation'] = stats[3]
    meta['ISO_RTC_transmittedPolarization'] = 'vertical'
    meta['ISO_RTC_receivedPolarization'] = 'vertical'
  raster = None

  raster = gdal.Open(fileName['original dem file'])
  band = raster.GetRasterBand(1)
  stats = band.GetStatistics(False, True)
  meta['ISO_RTC_DEM_minValue'] = stats[0]
  meta['ISO_RTC_DEM_maxValue'] = stats[1]
  meta['ISO_RTC_DEM_meanValue'] = stats[2]
  meta['ISO_RTC_DEM_standardDeviation'] = stats[3]
  raster = None
  meta['ISO_RTC_DEM_type'] = demSource

  raster = gdal.Open(fileName['layover shadow mask'])
  pixelCount = raster.RasterXSize * raster.RasterYSize
  raster = None
  lines = [line.rstrip() for line in open(fileName['layover shadow stats'])]
  for line in lines:
    if '0-  7:' in line:
      h = line.split(':')[1].split()
    if '8- 15:' in line:
      h += line.split(':')[1].split()
    if '16- 23:' in line:
      h += line.split(':')[1].split()
    if '24- 31:' in line:
      h += line.split(':')[1].split()
  noLayoverShadow = 0.0
  trueLayover = 0.0
  layover = 0.0
  trueShadow = 0.0
  shadow = 0.0
  pixelCount -= float(h[0])
  for ii in range(len(h)):
    if ii & 1:
      noLayoverShadow += float(h[ii]) / pixelCount
    if ii & 2:
      trueLayover += float(h[ii]) / pixelCount
    if ii & 4:
      layover += float(h[ii]) / pixelCount
    if ii & 8:
      trueShadow += float(h[ii]) / pixelCount
    if ii & 16:
      shadow += float(h[ii]) / pixelCount
  meta['ISO_RTC_noLayoverShadowPercentage'] = noLayoverShadow * 100.0
  meta['ISO_RTC_trueLayoverPercentage'] = trueLayover * 100.0
  meta['ISO_RTC_layoverPercentage'] = layover * 100.0
  meta['ISO_RTC_trueShadowPercentage'] = trueShadow * 100.0
  meta['ISO_RTC_shadowPercentage'] = shadow * 100.0

  raster = gdal.Open(fileName['incidence angle file'])
  band = raster.GetRasterBand(1)
  stats = band.GetStatistics(False, True)
  meta['ISO_RTC_incidenceMinValue'] = stats[0]
  meta['ISO_RTC_incidenceMaxValue'] = stats[1]
  meta['ISO_RTC_incidenceMeanValue'] = stats[2]
  meta['ISO_RTC_incidenceStandardDeviation'] = stats[3]
  raster = None

  # Processing
  lines = [line.rstrip() for line in open(fileName['mk_geo_radcal_0 log'])]
  for line in lines:
    if 'processing start:' in line:
      processing_time = \
        datetime.strptime(line.split(': ')[1].rsplit(' ',1)[0],
        '%a %b %d %H:%M:%S %Y ')
  meta['ISO_RTC_productCreationTime'] = processing_time.isoformat() + 'Z'

  return meta


def generate_product_dictionary(listFile):

  ### Read information from processing list file
  config = configparser.ConfigParser(allow_no_value=True)
  config.read(listFile)
  productType = config.sections()[0]

  ### Product: GAMMA RTC
  if productType == 'GAMMA RTC':
    meta = gamma_rtc_metadata(config)

  return meta


def get_metadata_values(metaFile, params):

  productParams = []

  parser = et.XMLParser(remove_blank_text=True)
  tree = et.parse(metaFile, parser)
  root = tree.getroot()
  for elem in root.getiterator():
    if len(elem.attrib) > 0:
      key = elem.attrib.keys()
      item = elem.attrib.items()
      (key, value) = item[0]
      if '}' in key:
        del elem.attrib[key]
        key = key.split('}')[1]
        elem.set(key, value)
    if not hasattr(elem.tag, 'find'):
      continue
    ii = elem.tag.find('}')
    if ii >= 0:
      elem.tag = elem.tag[ii+1:]
  objectify.deannotate(root, cleanup_namespaces=True)

  for ii in range(len(params)):
    try:
      param = tree.xpath('/'+params[ii])[0].text
      productParams.append(getValue(param))
    except:
      (element, attribute) = params[ii].rsplit('/',1)
      if attribute != 'value':
        elementAttribute = ('/'+element+'/@'+attribute)
        param = tree.xpath(elementAttribute)[0]
        productParams.append(getValue(param))
      else:
        param = tree.xpath('/'+element)[0].text
        productParams.append(getValue(param))

  return productParams


def product_dictionary2values(values, prodDict):

  prodValues = []
  for ii in range(len(values)):
    if ('{# ISO' or '{# UMM') in str(values[ii]):
      start = values[ii].index('{# ')
      stop = values[ii].index(' #}')
      parameter = values[ii][start:stop+3]
      key = parameter[3:-3]
      values[ii] = \
        values[ii].replace(parameter, str(prodDict[key]))
      try:
        if isinstance(prodDict[key], float) == True:
          values[ii] = float(values[ii])
        elif isinstance(prodDict[key], bool) == True:
          values[ii] = bool(values[ii])
        elif isinstance(prodDict[key], int) == True:
          values[ii] = int(values[ii])
      except:
        pass

    prodValues.append(values[ii])

  return prodValues


def granule_iso2umm_values(isoParams, isoValues, iso2umm, ummValues):

  (paramCount, _) = iso2umm.shape
  for kk in range(len(ummValues)):
    for ii in range(paramCount):
      ummParamString = ('{# UMM_%s #}' % iso2umm.at[ii, 'UMM'])
      if ummParamString == ummValues[kk]:
        isoIndex = isoParams.index(iso2umm.at[ii, 'ISO'])
        ummValues[kk] = isoValues[isoIndex]
  
  return ummValues


def properties2values(template, properties):

  values = []
  for param in template:
    if param in properties.keys():
      try:
        if isinstance(properties[param], float):
          values.append(float(properties[param]))
        elif isinstance(properties[param], bool):
          values.append(bool(properties[param]))
        elif isinstance(properties[param], int):
          values.append(int(properties[param]))
        else:
          values.append(properties[param])
      except:
        pass
    else:
      values.append(param)

  return values