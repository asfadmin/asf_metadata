#!/usr/bin/python3

"""Import packages"""
import argparse
from argparse import RawTextHelpFormatter
import os
import sys
import datetime
import json
from copy import deepcopy
from jinja2 import Template
from metadata.asf_metadata import project_template2dict


def iterate_over_dictionary(meta_input):
    """Recursively go through dictionary"""

    meta_output = deepcopy(meta_input)
    for key, value in meta_input.items():
        if isinstance(value, dict):
            meta_output[key] = iterate_over_dictionary(value)
        elif isinstance(value, list):
            meta_output[key] = []
            list_count = len(value)
            for i in range(list_count):
                meta_output[key].append(iterate_over_dictionary(value[i]))
        else:
            try:
                meta_output[key] = int(value)
            except ValueError:
                try:
                    meta_output[key] = float(value)
                except ValueError:
                    if value.startswith('[') and value.endswith(']'):
                        value_list = value.strip('][').split(', ')
                        list_count = len(value_list)
                        for i in range(list_count):
                            value_list[i] = float(value_list[i])
                        meta_output[key] = value_list
                    else:
                        meta_output[key] = value

    return meta_output


def generate_product_metadata(proc_file, excel_file, meta_level, meta_status,
    meta_file):
    """Generate product metadata"""

    ### Read processing JSON file into dictionary
    print(f'\nReading processing JSON file: {os.path.basename(proc_file)} ...')
    with open(proc_file, encoding='utf-8') as json_file:
        data = json.load(json_file)
    today = datetime.datetime.utcnow()
    data['timeSeries']['generation'] = today.isoformat() + 'Z'
    data['timeSeries']['status'] = meta_status
    data['timeSeries']['level'] = meta_level
    data['timeSeries']['generationYear'] = today.year
    data['timeSeries']['startYear'] = data['timeSeries']['stack']['startDate'][:4]
    data['timeSeries']['stopYear'] = data['timeSeries']['stack']['stopDate'][:4]

    ### Read project Excel spreadsheet template
    print('\nReading project metadata template: '
        f'{os.path.basename(excel_file)} ...')
    project_dict = project_template2dict(excel_file,
        'Project Metadata Structure', meta_level, data)
    project_dict_string = json.dumps(project_dict, indent=2)

    ### Put values into template
    print('\nGenerating metadata structure ...')
    jinja_template = Template(project_dict_string)
    meta_structure = jinja_template.render(data)
    meta = json.loads(meta_structure)
    meta = iterate_over_dictionary(meta)
    print(json.dumps(meta, indent=2))

    ### Write metadata structure to JSON file
    print('Writing metadata structure to JSON file: '
        f'{os.path.basename(meta_file)} ...')
    with open(meta_file, 'w', encoding='utf-8') as file:
        json.dump(meta, file, indent=2)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='generate_iso_metadata.py',
        description='Generate product metadata file from Excel spreadsheet',
        formatter_class=RawTextHelpFormatter)
    parser.add_argument('proc_file',
        help='file path of product processing JSON file')
    parser.add_argument('excel_file', help='name of product template file')
    parser.add_argument('meta_file',
        help='name of the output metadata JSON file')
    parser.add_argument('-meta_level', default='comprehensive',
        help='level of metadata (defaut: comprehensive)\n' \
            'choices: ["comprehensive", "short"]')
    parser.add_argument('-meta_status', default='draft',
        help='status of metadata (default: draft)\n' \
            'choices: ["draft", "final"]')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()

    generate_product_metadata(args.proc_file, args.excel_file, args.meta_level,
        args.meta_status, args.meta_file)
