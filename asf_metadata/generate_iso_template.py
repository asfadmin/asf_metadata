#!/usr/bin/python3

"""Import packages"""
import argparse
from argparse import RawTextHelpFormatter
import sys
from asf_metadata.metadata.asf_metadata import iso_template2lists, add_dem_lists, \
  iso_xml_structure, meta_xml_file, iso_dictionary_structure, meta_json_file


# Wrapper for the command line version
def generate_iso_template(excel_file, dem_file, iso_base):
    """Generate ISO template"""

    ### Extract ISO template from Excel spreadsheet
    print('\nExtracting ISO metadata template from Excel spreadsheet ...')
    (iso_template, iso_params, iso_values) = \
        iso_template2lists(excel_file, 'ISO Metadata Structure')
    if dem_file:
        (dem_template, dem_params, dem_values) = \
        iso_template2lists(dem_file, 'ISO Metadata Structure')
        (iso_template, iso_params, iso_values) = add_dem_lists(iso_template, \
        iso_params, iso_values, dem_template, dem_params, dem_values)

    ### Write ISO template to files
    print('Writing ISO metadata template to files ...')
    iso_structure = iso_xml_structure(excel_file, iso_template, iso_params,
        iso_values, True)
    meta_xml_file(iso_structure, iso_base + '.xml')
    iso_structure = iso_dictionary_structure(excel_file, dem_file, \
        iso_template, iso_params, iso_values)
    meta_json_file(iso_structure, iso_base + '.json')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='generate_iso_template.py',
        description='Generate ISO template file from Excel spreadsheet',
        formatter_class=RawTextHelpFormatter)
    parser.add_argument('excel_file',
        help='name of the Excel template spreadsheet')
    parser.add_argument('-dem', default=None,
        help='name of DEM template spreadsheet')
    parser.add_argument('iso_base',
        help='basename of the ISO metadata template file')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()

    generate_iso_template(args.excel_file, args.dem, args.iso_base)
