#!/usr/bin/python3

import argparse
from argparse import RawTextHelpFormatter
import sys
from asf_metadata import iso_template2lists, generate_product_dictionary, \
  product_dictionary2values, iso_xml_structure, meta_xml_file, \
  iso_dictionary_structure, meta_json_file 
  

def generate_iso_metadata(listFile, excelFile, isoBase):

  ### Read ISO template
  (isoTemplate, isoParams, isoValues) = \
    iso_template2lists(excelFile, 'ISO Metadata Structure')

  ### Generate product dictionary
  isoProdDict = generate_product_dictionary(listFile)
  
  ### Map product dictionary to values
  isoProdValues = product_dictionary2values(isoValues, isoProdDict)

  ### Write ISO metadata to XML
  isoFile = isoBase + '.iso.xml'
  print('Writing ISO metadata structure to XML file ({0}) ...'.format(isoFile))
  isoStructure = \
    iso_xml_structure(excelFile, isoTemplate, isoParams, isoProdValues, True)
  meta_xml_file(isoStructure, isoFile)

  ### Write ISO metadata to JSON
  isoFile = isoBase + '.iso.json'
  print('Writing ISO metadata structure to JSON file ({0}) ...' \
    .format(isoFile))
  isoStructure = \
    iso_dictionary_structure(excelFile, isoTemplate, isoParams, isoProdValues)
  meta_json_file(isoStructure, isoFile)


if __name__ == '__main__':

  parser = argparse.ArgumentParser(prog='generate_iso_metadata.py',
    description='Generate ISO metadata file from Excel spreadsheet',
    formatter_class=RawTextHelpFormatter)
  parser.add_argument('listFile', help='name of processing list file')
  parser.add_argument('excelFile', 
    help='name of the Excel template spreadsheet')
  parser.add_argument('isoBase', help='basename of the ISO XML metadata file')
  if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
  args = parser.parse_args()

  generate_iso_metadata(args.listFile, args.excelFile, args.isoBase)