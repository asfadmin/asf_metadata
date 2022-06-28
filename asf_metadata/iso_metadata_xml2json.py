#!/usr/bin/python3

import argparse
from argparse import RawTextHelpFormatter
import sys
from asf_metadata.metadata.asf_metadata import iso_template2lists, add_dem_lists, \
  iso_dictionary_structure, get_metadata_values, meta_json_file


def iso_metadata_xml2json(excelFile, demFile, xmlFile, jsonFile):

  ### Read ISO template
  (isoTemplate, isoParams, isoValues) = \
    iso_template2lists(excelFile, 'ISO Metadata Structure')
  if demFile:
    (demTemplate, demParams, demValues) = \
      iso_template2lists(demFile, 'ISO Metadata Structure')
    (isoTemplate, isoParams, isoValues) = add_dem_lists(isoTemplate, \
      isoParams, isoValues, demTemplate, demParams, demValues)

  ### Extract values from ISO XML
  isoProdValues = get_metadata_values(xmlFile, isoParams)

  ### Write ISO metadata to JSON
  isoStructure = iso_dictionary_structure(excelFile, demFile, isoTemplate, \
    isoParams, isoProdValues)
  meta_json_file(isoStructure, jsonFile)


if __name__ == '__main__':

  parser = argparse.ArgumentParser(prog='iso_metadata_xml2json.py',
    description='Convert ISO metadata from XML to JSON format',
    formatter_class=RawTextHelpFormatter)
  parser.add_argument('excelFile', help='name of the Excel template spreadsheet')
  parser.add_argument('-dem', default=None,
    help='name of DEM template spreadsheet')
  parser.add_argument('xmlFile', help='name of the ISO metadata XML file')
  parser.add_argument('jsonFile', help='name of the ISO metadata JSON file')
  if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
  args = parser.parse_args()

  iso_metadata_xml2json(args.excelFile, args.dem, args.xmlFile, args.jsonFile)
