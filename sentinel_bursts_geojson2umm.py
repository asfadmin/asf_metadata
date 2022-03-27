#!/usr/bin/python3

import argparse
from argparse import RawTextHelpFormatter
import os
import sys
import json
from metadata.asf_metadata import iso_template2lists, properties2values, \
  umm_dictionary_structure, meta_json_file


def sentinel_bursts_geojson2umm(excelFile, jsonFile, outDir):

  ### Check output directory
  if not os.path.exists(outDir):
    os.mkdir(outDir)

  ### Get values for template placeholders
  with open(jsonFile) as inF:
    data = json.load(inF)

  ### Loop through all bursts from the granule
  features = data['features']
  for feature in features:

    ### Read UMM-G template
    (ummTemplate, ummParams, ummValues) = \
      iso_template2lists(excelFile, 'UMM-G Structure')
    
    ### Extract UMM values
    properties = feature['properties']
    properties['insertDate'] = '{# UMM_insertDate #}'
    properties['hasLand'] = '{ UMM_hasLand #}'
    properties['crossesDateLine'] = '{# UMM_crossesDateLine #}'
    properties['landFraction'] = '{# UMM_landFraction #}'
    ummDict = {}
    for key in properties.keys():
      ummKey = ('{# UMM_%s #}' % key)
      ummDict[ummKey] = properties[key]
    properties = ummDict
    ummProdValues = properties2values(ummValues, properties)
    burstID = properties['{# UMM_absoluteBurstID #}']

    ### Write UMM-G to file
    ummFile = ('%s-burst-umm.json' % burstID)
    print('Writing UMM-G metadata to JSON file ({0}) ...'.format(ummFile))
    ummStructure = \
      umm_dictionary_structure(ummTemplate, ummParams, ummProdValues)
    meta_json_file(ummStructure, os.path.join(outDir, ummFile))


if __name__ == '__main__':

  parser = argparse.ArgumentParser(prog='sentinel_bursts_geojson2umm.py',
    description='Convert Sentinel-1 bursts metadata from GeoJSON to ' \
    'UMM JSON format', formatter_class=RawTextHelpFormatter)
  parser.add_argument('excelFile', metavar='<Excel template>',
    help='name of the Excel template spreadsheet')
  parser.add_argument('jsonFile', metavar='<GeoJSON file>',
    help='name of the Sentinel-1 bursts metadata GeoJSON file')
  parser.add_argument('outDir', metavar='<output directory>',
    help='name of the output directory')
  if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
  args = parser.parse_args()

  sentinel_bursts_geojson2umm(args.excelFile, args.jsonFile, args.outDir)