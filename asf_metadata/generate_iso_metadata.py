#!/usr/bin/python3

"""Import packages"""
import argparse
from argparse import RawTextHelpFormatter
import os
import sys
from asf_metadata.metadata.asf_metadata import iso_template2lists, add_dem_lists, \
  generate_product_dictionary, product_dictionary2values, iso_xml_structure, \
  meta_xml_file, iso_dictionary_structure, meta_json_file, \
  clean_json_structure, clean_xml_structure


def generate_iso_metadata(meta_file, log_file, iso_base, product_type,
    data_type, dem_type):
    """Generate ISO metadata"""

    ### Work out the location of the templates and their files
    template_path = \
        os.path.join(os.path.dirname(os.path.realpath(__file__)), 'templates')
    data_template_file = f"{data_type}_{product_type}_metadata.xlsx"
    excel_file = os.path.join(template_path, data_type, data_template_file)
    print(f'\nData template: {excel_file}')
    if dem_type:
        dem_template_file = f"dem_{dem_type}_metadata.xlsx"
        dem_file = os.path.join(template_path, 'dem', dem_template_file)
        print(f'DEM template: {dem_file}')
    else:
        dem_file = None

    ### Read ISO template
    (iso_template, iso_params, iso_values) = \
        iso_template2lists(excel_file, 'ISO Metadata Structure')
    if dem_type:
        (dem_template, dem_params, dem_values) = \
            iso_template2lists(dem_file, 'ISO Metadata Structure')
        (iso_template, iso_params, iso_values) = add_dem_lists(iso_template, \
            iso_params, iso_values, dem_template, dem_params, dem_values)
    else:
        dem_params = None
        dem_template = None

    ### Generate product dictionary
    iso_prod_dict = generate_product_dictionary(product_type,
        data_type, meta_file, log_file)

    ### Map product dictionary to values
    iso_prod_values = product_dictionary2values(iso_values, iso_prod_dict)

    ### Write ISO metadata to XML
    iso_file = iso_base + '.iso.xml'
    print(f'Writing ISO metadata structure to XML file ({iso_file}) ...')
    iso_structure = iso_xml_structure(excel_file, iso_template, iso_params, \
        iso_prod_values, True)
    iso_structure = clean_xml_structure(iso_file)
    meta_xml_file(iso_structure, iso_file)

    ### Write ISO metadata to JSON
    iso_file = iso_base + '.iso.json'
    print(f'Writing ISO metadata structure to JSON file ({iso_file}) ...')
    iso_structure = iso_dictionary_structure(excel_file, dem_file, \
        iso_template, iso_params, iso_prod_values)
    iso_structure = clean_json_structure(iso_structure)
    meta_json_file(iso_structure, iso_file)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='generate_iso_metadata.py',
        description='Generate ISO metadata file from Excel spreadsheet',
        formatter_class=RawTextHelpFormatter)
    parser.add_argument('meta_file', help='file path of product metadata ' \
        'such as a "manifest.safe" or leader file')
    parser.add_argument('log_file', help='name of RTC processing log file')
    parser.add_argument('iso_base', help='basename of the ISO XML metadata file')
    parser.add_argument('-product', default='gamma_rtc',
        help='name of the product type (default: gamma_rtc)\n' \
            'choices: ["gamma_rtc"]')
    parser.add_argument('-data', default='sentinel',
        help='name of the data source (default: sentinel)\n' \
            'choices: ["sentinel"]')
    parser.add_argument('-dem', default=None,
        help='DEM type used for terrain correction\n' \
            'choices: ["copernicus"]')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()

    generate_iso_metadata(args.meta_file, args.log_file, args.iso_base,
        args.product, args.data, args.dem)
