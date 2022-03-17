#!/usr/bin/python3

import argparse
from argparse import RawTextHelpFormatter
import sys
from metadata.asf_metadata import iso_template2lists, add_dem_lists, \
  generate_product_dictionary, product_dictionary2values, iso_xml_structure, \
  meta_xml_file, iso_dictionary_structure, meta_json_file, \
  cleanJSONstructure, cleanXMLstructure
  

def generate_iso_metadata(productType, logFile, excelFile, demFile, isoBase):

  if demFile is not None:
    print('DEM template: {0}'.format(demFile))

  ### Read ISO template
  (isoTemplate, isoParams, isoValues) = \
    iso_template2lists(excelFile, 'ISO Metadata Structure')
  if demFile:
    (demTemplate, demParams, demValues) = \
      iso_template2lists(demFile, 'ISO Metadata Structure')
    (isoTemplate, isoParams, isoValues) = add_dem_lists(isoTemplate, \
      isoParams, isoValues, demTemplate, demParams, demValues)
  else:
    demParams = None
    demTemplate = None

  ### Generate product dictionary
  isoProdDict = generate_product_dictionary(productType, logFile)
  
  ### Map product dictionary to values
  isoProdValues = product_dictionary2values(isoValues, isoProdDict)

  ### Write ISO metadata to XML
  isoFile = isoBase + '.iso.xml'
  print('Writing ISO metadata structure to XML file ({0}) ...'.format(isoFile))
  isoStructure = \
    iso_xml_structure(excelFile, isoTemplate, isoParams, isoProdValues, True)
  meta_xml_file(isoStructure, isoFile)
  isoStructure = cleanXMLstructure(isoFile)
  meta_xml_file(isoStructure, isoFile)

  ### Write ISO metadata to JSON
  isoFile = isoBase + '.iso.json'
  print('Writing ISO metadata structure to JSON file ({0}) ...' \
    .format(isoFile))
  isoStructure = iso_dictionary_structure(excelFile, demFile, isoTemplate, \
    isoParams, isoProdValues)
  isoStructure = cleanJSONstructure(isoStructure)
  meta_json_file(isoStructure, isoFile)


if __name__ == '__main__':

  parser = argparse.ArgumentParser(prog='generate_iso_metadata.py',
    description='Generate ISO metadata file from Excel spreadsheet',
    formatter_class=RawTextHelpFormatter)
  parser.add_argument('productType', help='name of the product type')
  parser.add_argument('logFile', help='name of processing log file')
  parser.add_argument('excelFile', 
    help='name of the Excel template spreadsheet')
  parser.add_argument('-dem', default=None,
    help='name of DEM template spreadsheet')
  parser.add_argument('isoBase', help='basename of the ISO XML metadata file')
  if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
  args = parser.parse_args()

  generate_iso_metadata(args.productType.upper(), args.logFile, args.excelFile,
    args.dem, args.isoBase)
