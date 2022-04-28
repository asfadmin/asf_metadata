#!/usr/bin/python3

"""Import packages"""
import argparse
from argparse import RawTextHelpFormatter
import os
import sys
import json
from metadata.asf_metadata import iso_template2lists, properties2values, \
  umm_dictionary_structure, meta_json_file


def sentinel_bursts_geojson2umm(excel_file, json_file, out_dir):
    """Convert burst metadata from GeoJSON to UMM-G format"""

    ### Check output directory
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    ### Get values for template placeholders
    with open(json_file, encoding='utf-8') as file:
        data = json.load(file)

    ### Loop through all bursts from the granule
    features = data['features']
    for feature in features:

        ### Read UMM-G template
        (umm_template, umm_params, umm_values) = \
            iso_template2lists(excel_file, 'UMM-G Structure')

        ### Extract UMM values
        properties = feature['properties']
        properties['insertDate'] = '{# UMM_insertDate #}'
        properties['hasLand'] = '{# UMM_hasLand #}'
        properties['crossesDateLine'] = '{# UMM_crossesDateLine #}'
        properties['landFraction'] = '{# UMM_landFraction #}'
        umm_dict = {}
        for key in properties.keys():
            umm_key = ('{# UMM_%s #}' % key)
            umm_dict[umm_key] = properties[key]
        properties = umm_dict
        umm_prod_values = properties2values(umm_values, properties)
        burst_id = properties['{# UMM_absoluteBurstID #}']

        ### Write UMM-G to file
        umm_file = f'{burst_id}-burst-umm.json'
        print(f'Writing UMM-G metadata to JSON file ({umm_file}) ...')
        umm_structure = \
            umm_dictionary_structure(umm_template, umm_params, umm_prod_values)
        meta_json_file(umm_structure, os.path.join(out_dir, umm_file))


if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='sentinel_bursts_geojson2umm.py',
        description='Convert Sentinel-1 bursts metadata from GeoJSON to ' \
        'UMM JSON format', formatter_class=RawTextHelpFormatter)
    parser.add_argument('excel_file', metavar='<Excel template>',
        help='name of the Excel template spreadsheet')
    parser.add_argument('json_file', metavar='<GeoJSON file>',
        help='name of the Sentinel-1 bursts metadata GeoJSON file')
    parser.add_argument('out_dir', metavar='<output directory>',
        help='name of the output directory')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()

    sentinel_bursts_geojson2umm(args.excel_file, args.json_file, args.out_dir)
