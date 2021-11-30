#!/usr/bin/python3

import argparse
from argparse import RawTextHelpFormatter
import sys
from asf_metadata import iso_template2lists, meta_json_file, get_metadata_values


def iso_metadata_xml2json(excelFile, xmlFile, jsonFile):

  ### Read ISO template
  (_, _, isoParams, _) = iso_template2lists(excelFile)

  ### Extract values from ISO XML
  productValues = get_metadata_values(xmlFile, isoParams)

  ### Generate ISO metadata to JSON
  meta_json_file(isoParams, productValues, jsonFile)


if __name__ == '__main__':

  parser = argparse.ArgumentParser(prog='iso_metadata_xml2json.py',
    description='Convert ISO metadata from XML to JSON format',
    formatter_class=RawTextHelpFormatter)
  parser.add_argument('excelFile', help='name of the Excel template spreadsheet')
  parser.add_argument('xmlFile', help='name of the ISO metadata XML file')
  parser.add_argument('jsonFile', help='name of the ISO metadata JSON file')
  if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
  args = parser.parse_args()

  iso_metadata_xml2json(args.excelFile, args.xmlFile, args.jsonFile)
