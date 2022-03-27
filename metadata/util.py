#!/usr/bin/python3

from itertools import groupby
from operator import itemgetter
import pandas as pd
import configparser
from sqlalchemy import create_engine


def splitElement(k, v, out):

  k, *rest = k.split('/', 1)
  if rest:
    splitElement(rest[0], v, out.setdefault(k, {}))
  else:
    out[k] = v


def nestDict(flatDict):

  result = {}
  for k, v in flatDict.items():
    splitElement(k, v, result)

  return result


def isFloat(value):

  try:
    float(value)
    return True
  except ValueError:
    return False


def isInt(value):

  try:
    int(value)
    return True
  except ValueError:
    return False


def getValue(valueStr):

  if isInt(valueStr) == True:
    outValue = int(valueStr)
  elif isFloat(valueStr) == True:
    outValue = float(valueStr)
  else:
    outValue = valueStr.strip()
  return outValue


def rreplace(string, find, replace):
    reversed = string[::-1]
    replaced = reversed.replace(find[::-1], replace[::-1], 1)
    return replaced[::-1]


def uniqueList(strList):
	
  outList = []
  for x in strList:
    if x not in outList:
      outList.append(x)
  return outList


def meta2dict(metaFile):
	
  ### Extract parameters and values from metadata file
  metaParams = []
  metaValues = []
  ii = 0
  section = ''
  lines = [line.rstrip() for line in open(metaFile)]
  for line in lines:
    if '{' in line:
      group = line.split(' {')[0].lstrip()
      section = group
      if section == 'vector':
        ii += 1
    elif ':' in line:
      param = line.split(':')[0].lstrip()
      if section == 'vector':
        metaParams.append('state/vector[{0}]/{1}'.format(ii, param))
      elif section == 'param':
        metaParams.append('projection/param/{0}'.format(param))
      elif section == 'utm':
        metaParams.append('projection/param/utm/{0}'.format(param))
      elif section == 'band_stats':
        metaParams.append('statistics/band_stats/{0}'.format(param))
      else:
        if param == 'meta_version':
          metaParams.append(param)
        else:
          metaParams.append(section + '/' + param)
      value = line.split(': ')[1].split('#')[0].rstrip()
      if isInt(value) == True:
        metaValues.append(int(value))
      elif isFloat(value) == True:
        metaValues.append(float(value))
      else:
        metaValues.append(value)
    elif '}' in line:
      section = ''
    else:
      continue

  ### Establish flattened dictionary
  flatDict = dict(zip(metaParams, metaValues))

  ### Convert flattened into nested dictionary
  metaDict = nestDict(flatDict)

  return metaDict


def mergeDict(a, b):

  key = None
  if a is None or isinstance(a, (str, int, float)):
    a = b
  elif isinstance(a, list):
    if isinstance(b, list):
      #print('a: {0}'.format(a))
      #print('b: {0}'.format(b))
      if len(a[0]) == len(b[0]):
        #a.append(b[0])
        a[0] = mergeDict(a[0], b[0])
      elif next(iter(b[0])) not in a[len(a)-1].keys():
        #print('elif')
        a[len(a)-1] = mergeDict(a[len(a)-1], b[0])
      else:
        #print('else')
        a.extend(b)
    else:
      a.append(b)
  elif isinstance(a, dict):
    if isinstance(b, dict):
      for key in b:
        if key in a:
          a[key] = mergeDict(a[key], b[key])
        else:
          a[key] = b[key]
  
  return a


def getConsecutiveList(metaList):

  for _, g in groupby(enumerate(metaList), lambda x: x[0]-x[1]):
    yield list(map(itemgetter(1), g))


def getLevelParamsList(metaTemplate, metaAttributes, level):

  levelCount = 0
  rowCount = len(metaTemplate)
  for ii in range(rowCount):
    metaElement = metaTemplate[ii].split('/')
    levelCount = \
      len(metaElement) if len(metaElement) > levelCount else levelCount

  df = pd.DataFrame(index=range(rowCount), columns=range(levelCount))
  for ii in range(rowCount):
    metaElement = metaTemplate[ii].split('/')
    for kk in range(len(metaElement)):
      df.at[ii,kk] = metaElement[kk]
    if metaAttributes:
      if isinstance(metaAttributes[ii], str) and level == kk+1:
        df.at[ii,kk+1] = metaAttributes[ii]
  rowIndexList = list(df.count(axis=1) - 1)
  indices = [i for i, x in enumerate(rowIndexList) if x == level]
  indices = list(df.filter(items=indices, axis=0).index)
  levelParams = []
  for index in indices:
    param = '/'.join(list(df.loc[index].dropna()))
    if param.endswith(']'):
      levelParams.append(param.rsplit('[',1)[0])
    else:
      levelParams.append(param)
  
  return levelParams


def getParamsDataframe(metaParams):

  levelCount = 0
  rowCount = len(metaParams)
  for ii in range(rowCount):
    metaElement = metaParams[ii].split('/')
    levelCount = \
      len(metaElement) if len(metaElement) > levelCount else levelCount

  df = pd.DataFrame(index=range(rowCount), columns=range(levelCount-1),
    dtype=int)
  for ii in range(rowCount):
    metaElement = metaParams[ii].split('/')
    for kk in range(levelCount):
      if kk < len(metaElement):
        if metaElement[kk].endswith(']'):
          df.at[ii,kk] = int(metaElement[kk].split('[')[1][:-1])
        else:
          df.at[ii,kk] = 0
      else:
        df.at[ii,kk] = None

  return df


def getDictListParams(dictList, oldLevel, newLevel):

  dictParents = []
  dictParams = []
  for dictElement in dictList:
    dictParents.append(dictElement['path'].rsplit('/',oldLevel-newLevel)[0])
    dictParams.append(dictElement['path'])

  return (dictParents, dictParams)


def params2dictList(df, metaParams, metaValues, level):

  dictList = []
  rowIndex = list(df.count(axis=1) - 1)
  colIndex = list(df.sum(axis=0))
  for ii in range(len(colIndex)):
    if colIndex[ii] > 0:
      rowIndex.append(ii)
  rowIndex = pd.Series(rowIndex)
  dfColumn = df[level].dropna()
  rowLevelIndex = rowIndex[rowIndex==level]
  dfColumn = dfColumn.filter(items=list(rowLevelIndex.index))
  dfColumnLength = len(dfColumn)
  if dfColumnLength > 0:
    minList = int(dfColumn.min())
    maxList = int(dfColumn.max())
    for listIndex in range(minList, maxList+1):
      dfList = dfColumn[dfColumn == listIndex]
      dfIndexList = list(getConsecutiveList(dfList.index.tolist()))
      for dfIndices in dfIndexList:
        for dfIndex in dfIndices:
          metaElement = metaParams[dfIndex].split('/')
          path = metaParams[dfIndex]
          parent = {}
          parent['level'] = len(path.split('/'))-1
          parent['path'] = path
          parent['type'] = 'value'
          parent['key'] = metaElement[level]
          parent['list'] = listIndex
          parent['dict'] = metaValues[dfIndex]
          dictList.append(parent)

  return dictList


def upgradeDictionary2level(oldDictList, oldLevel, newLevel):

  newDictList = []
  (dictParents, _) = \
    getDictListParams(oldDictList, oldLevel, newLevel)
  parents = list(dict.fromkeys(dictParents))
  
  for parent in parents:
    indices = [i for i, x in enumerate(dictParents) if x == parent]
    for index in indices:
      parentDict = {}
      levelDiff = oldDictList[index]['level'] - newLevel
      path = oldDictList[index]['path']
      oldDict = oldDictList[index]['dict']
      for _ in range(levelDiff):
        key = path.rsplit('/',1)[1]
        if ('Code' in key) and (key != 'postalCode'):
          value = oldDict
          oldDict = {}
          if 'EOS' in key:
            codeURL = 'http://earthdata.nasa.gov/metadata/resources/' \
              'Codelist.xml#' + key
          else:
            codeURL = 'https://cdn.earthdata.nasa.gov/iso/resources/' \
              'Codelist/gmxCodelists.xml#' + key
          oldDict['@codeList'] = codeURL
          oldDict['@codeListValue'] = value
          oldDict['value'] = value
        newDict = {}
        newDict[key.split('[')[0]] = oldDict
        path = path.rsplit('/',1)[0]
        oldDict = newDict
      parentDict['level'] = newLevel
      if path.endswith(']'):
        listIndex = int(path.rsplit('[',1)[1][:-1])
        path = path.rsplit('[',1)[0]
        parentDict['path'] = path
        parentDict['key'] = path.rsplit('/',1)[1]
        parentDict['list'] = listIndex
        parentDict['type'] = 'list'
      else:
        if isinstance(newDict, dict):
          parentDict['path'] = path
        elif isinstance(newDict, (int, float, str)):
          parentDict['path'] = path.rsplit('/',1)[0]
        if '/' in path:
          parentDict['key'] = path.rsplit('/',1)[1]
        else:
          parentDict['key'] = path
        parentDict['list'] = 0
        parentDict['type'] = 'single'
        parentDict['dict'] = newDict
      if levelDiff > 0:
        parentDict['dict'] = newDict
      else:
        parentDict['dict'] = oldDictList[index]['dict']
      newDictList.append(parentDict)

  return newDictList


def mergeDictionaryParams(oldDictList, levelParams, oldLevel, newLevel):

  mergedDictList = []
  (_, dictParams) = \
    getDictListParams(oldDictList, oldLevel, newLevel)
  params = list(dict.fromkeys(dictParams))
  params = [i for i, g in groupby(levelParams)]

  for param in params:
    indices = [i for i, x in enumerate(dictParams) if x == param]
    listIndices = []
    for index in indices:
      listIndices.append(oldDictList[index]['list'])
    if len(listIndices) > 0:
      maxIndex = max(listIndices)
    else:
      maxIndex = 0

    ### list case
    if maxIndex > 0:
      dictIndices = [[] for i in range(maxIndex)]
      for index in indices:
        dictIndices[oldDictList[index]['list']-1].append(index)
      numIndices = len(dictIndices[0])
      combineDict = []
      for dictIndex in dictIndices:
        oldDict = oldDictList[dictIndex[0]]
        if numIndices > 1:
          listDict = {}
          for listIndex in dictIndex:
            oldDict = oldDictList[listIndex]
            key, value = next(iter(oldDict['dict'].items()))
            listDict[key] = value
          combineDict.append(listDict)
        else:
          combineDict.append(oldDict['dict'])

    ### single value case
    else:
      combineDict = {}
      dictType = None
      for index in indices:
        oldDict = oldDictList[index]
        if isinstance(oldDict['dict'], dict):
          key, value = next(iter(oldDict['dict'].items()))
          combineDict[key] = value
          dictType = 'dict'
        elif isinstance(oldDict['dict'], (float, int, str)):
          if dictType == 'dict':
            combineDict['value'] = oldDict['dict']
          else:
            combineDict = oldDict['dict']

    mergedDict = {}
    mergedDict['level'] = oldDict['level']
    if oldDict['key'] == 'type' and \
      'acquisitionInformation' in oldDict['path']:
      mergedDict['path'] = param
      mergedDict['key'] = oldDict['key']
    elif oldDict['key'] in ['id','href','nilReason','type','uom'] and \
      'additionalAttribute' not in oldDict['path']:
      mergedDict['path'] = rreplace(param, '/','/@')
      mergedDict['key'] = '@' + oldDict['key']
    else:
      mergedDict['path'] = param
      mergedDict['key'] = oldDict['key']
    mergedDict['list'] = oldDict['list']
    mergedDict['type'] = 'single'
    mergedDict['dict'] = combineDict
    mergedDictList.append(mergedDict)

  return mergedDictList


def create_postgis_engine(configFile):

  config = configparser.ConfigParser()
  config.read(configFile)
  connection_string = \
    "postgresql://" + config.get('postgres', 'user') + \
    ":'" + config.get('postgres', 'pass') + "'" + \
    "@" + config.get('postgres', 'host') + ":5432" + \
    "/" + config.get('postgres', 'db')
  engine = create_engine(connection_string)

  return engine
