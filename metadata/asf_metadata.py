#!/usr/bin/python3

import lxml.etree as et
from lxml import objectify
import json
import pandas as pd
import numpy as np
from datetime import datetime
import collections
import re
import os
from osgeo import gdal, osr
import glob
from metadata.util import getValue, rreplace, uniqueList, \
  getParamsDataframe, params2dictList, upgradeDictionary2level, \
  mergeDictionaryParams, getLevelParamsList
from metadata.sentinelMetadata import parseLine, get_latlon_extent, \
  readManifestFile, readAnnotationFile


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


def update_template(template, listParams):

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
    if len(listParams) > 0:
      isoTemplate = update_template(isoTemplate, listParams)

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
    json.dump(metaStructure, outF, indent=2)


def meta_xml_file(metaStructure, xmlFile):

  ### Write XML file
  with open(xmlFile, 'wb') as outF:
    outF.write(et.tostring(metaStructure, xml_declaration=True, encoding='utf-8',
      pretty_print=True))


def generate_product_dictionary(productType, dataSource, metaFile, logFile):

  if productType == 'GAMMA RTC':
    meta = gammaRTClog2meta(dataSource, metaFile, logFile)

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

  dictKeys = list(prodDict.keys())
  prodValues = []
  for ii in range(len(values)):
    keys = []
    for kk in range(len(dictKeys)):
      if str(values[ii]).count(dictKeys[kk]) > 0:
        keys.append(dictKeys[kk])
    for key in keys:
      start = values[ii].index('{# ')
      stop = values[ii].index(' #}')
      parameter = values[ii][start:stop+3]
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


def cleanJSONstructure(isoStructure):

  meta = \
    isoStructure['DS_Series']['composedOf']['DS_DataSet']['has']['MI_Metadata']

  identificationInfo = meta['identificationInfo']
  newIdentification = []
  for ii in range(len(identificationInfo)):
    value = identificationInfo[ii]['MD_DataIdentification']['citation'] \
      ['CI_Citation']['title']['FileName']
    if not value.startswith('{# ISO'):
      newIdentification.append(identificationInfo[ii])
  meta['identificationInfo'] = newIdentification

  contentInfo = meta['contentInfo']
  newContent = []
  for ii in range(len(contentInfo)):
    miBand = contentInfo[ii]['MD_CoverageDescription']['dimension']['MI_Band']
    if 'maxValue' in miBand:
      value = miBand['maxValue']['Real']
    elif 'otherPropertyType' in miBand:
      value = miBand['otherProperty']['Record']['additionalAttributes'] \
        ['EOS_AdditionalAttributes']['additionalAttribute'][0] \
        ['EOS_AdditionalAttribute']['value']['CharacterString']
    if type(value) == float:
      newContent.append(contentInfo[ii])
  meta['contentInfo'] = newContent

  return isoStructure


def cleanXMLstructure(isoFile):

  ns_gmd = {'gmd' : 'http://www.isotc211.org/2005/gmd'}
  ns_gco = {'gco' : 'http://www.isotc211.org/2005/gco'}
  ns_eos = {'eos' : 'http://earthdata.nasa.gov/schema/eos'}
  ns_gmi = {'gmi' : 'http://www.isotc211.org/2005/gmi'}
  ns_gmx = {'gmx' : 'http://www.isotc211.org/2005/gmx'}
  ns = dict(
    list(ns_gmd.items()) +
    list(ns_gco.items()) +
    list(ns_eos.items()) +
    list(ns_gmi.items()) +
    list(ns_gmx.items())
  )

  parser = et.XMLParser(remove_blank_text=True)
  doc = et.parse(isoFile, parser)
  meta = doc.xpath('/gmd:DS_Series/gmd:composedOf/gmd:DS_DataSet' \
    '/gmd:has/gmi:MI_Metadata', namespaces=ns)[0]
  identificationInfoCount = 0
  contentInfoCount = 0
  remove = []
  for ii in range(len(meta)):
    if 'identificationInfo' in meta[ii].tag:
      fileName = doc.xpath('/gmd:DS_Series/gmd:composedOf/gmd:DS_DataSet' \
        '/gmd:has/gmi:MI_Metadata/gmd:identificationInfo[{0}]' \
        '/gmd:MD_DataIdentification/gmd:citation/gmd:CI_Citation/gmd:title' \
        '/gmx:FileName'.format(identificationInfoCount+1), 
        namespaces=ns)[0].text
      if fileName.startswith('{# ISO'):
        remove.append(ii)
      identificationInfoCount += 1
    elif 'contentInfo' in meta[ii].tag:
      miBand = doc.xpath('/gmd:DS_Series/gmd:composedOf/gmd:DS_DataSet' \
        '/gmd:has/gmi:MI_Metadata/gmd:contentInfo[{0}]' \
        '/gmd:MD_CoverageDescription/gmd:dimension/gmi:MI_Band' \
        .format(contentInfoCount+1), namespaces=ns)[0]
      if 'maxValue' in miBand[0].tag:
        value = doc.xpath('/gmd:DS_Series/gmd:composedOf/gmd:DS_DataSet' \
        '/gmd:has/gmi:MI_Metadata/gmd:contentInfo[{0}]' \
        '/gmd:MD_CoverageDescription/gmd:dimension/gmi:MI_Band/gmd:maxValue' \
        '/gco:Real'.format(contentInfoCount+1), namespaces=ns)[0].text
      elif 'otherPropertyType' in miBand[0].tag:
        value = doc.xpath('/gmd:DS_Series/gmd:composedOf/gmd:DS_DataSet' \
        '/gmd:has/gmi:MI_Metadata/gmd:contentInfo[{0}]' \
        '/gmd:MD_CoverageDescription/gmd:dimension/gmi:MI_Band' \
        '/eos:otherProperty/gco:Record/eos:additionalAttributes' \
        '/eos:EOS_AdditionalAttributes/eos:additionalAttribute[1]' \
        '/eos:EOS_AdditionalAttribute/eos:value/eos:CharacterString' \
        .format(contentInfoCount+1), namespaces=ns)[0].text
      try:
        float(value)
      except ValueError:
        remove.append(ii)
      contentInfoCount += 1
 
  for ii in range(len(remove)-1,-1,-1):
    del meta[remove[ii]]
  
  return doc


def gammaRTClog2meta(dataSource, metaFile, logFile):

  ns_safe = {'safe': 'http://www.esa.int/safe/sentinel-1.0'}
  ns = dict(
    list(ns_safe.items())
  )

  ### Determine parent directory of log file
  parentDir = os.path.dirname(os.path.abspath(os.path.basename(metaFile)))

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

  if dataSource == 'SENTINEL':

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

    ### Examine product directory for various files and calculate stats if exists
    polarization = meta['ISO_RTC_inputGranule'][14:16]
    if polarization == 'SH' or polarization == 'DH':
      mainPol = 'HH'
      polString = 'horizontal'
    elif polarization == 'SV' or polarization == 'DV':
      mainPol = 'VV'
      polString = 'vertical'

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
