#!/usr/bin/python3

"""Import packages"""
import json
from datetime import datetime
import collections
import re
import os
import glob
import lxml.etree as et
from lxml import objectify
import pandas as pd
import numpy as np
from osgeo import gdal, osr
import flatdict
from metadata.util import get_value, rreplace, unique_list, \
    get_params_dataframe, params2dict_list, upgrade_dictionary2level, \
    merge_dictionary_params, get_level_params_list, dataframe2template
from metadata.sentinel_metadata import parse_line, get_latlon_extent, \
    read_manifest_file, read_annotation_file


def meta_project2epsg(meta_dict):
    """Determine EPSG code from projection type"""

    projection = meta_dict['projection']['type']
    if projection == 'UNIVERSAL_TRANSVERSE_MERCATOR':
        hemisphere = meta_dict['projection']['hem']
        zone = int(meta_dict['projection']['param']['utm']['zone'])
        if hemisphere == 'N':
            epsg = 32600 + zone
        else:
            epsg = 32700 + zone
    return epsg


def update_template(template, list_params):
    """Update template"""

    count = len(template)
    for list_param in list_params:
        indices = [i for i, x in enumerate(template) if x == list_param]
        indices.append(count)
        for i in range(len(indices)-1):
            for k in range(indices[i], indices[i+1]):
                if list_param in template[k]:
                    new_element = f'{list_param}[{i+1}]'
                    template[k] = \
                        rreplace(template[k], list_param, new_element)

    return template


def get_duplicate_indices(template, list_param):
    """Get start and stop indices of duplicates in the template"""

    param_list = []
    param_list.append(list_param)
    duplicates = [substring for substring in template
        if all((element in substring) for element in param_list)]
    duplicate_count = len(duplicates)
    indices = []
    for i in range(duplicate_count):
        indices.append(template.index(duplicates[i]))

    return indices


def add_duplicates2metadata(meta, list_params, list_counts):
    """Add duplicates to metadata"""

    ### Check dimensions
    (row_count, column_count) = meta.shape
    level_count = column_count - 6

    ### Convert dataframe to template
    template = dataframe2template(meta, level_count)

    ### Add duplicate parameters
    new_meta = pd.DataFrame()
    duplicate_count = len(list_params)
    for i in range(row_count):
        if '[i]' not in str(meta.at[i,'Value']) and \
            meta.at[i,'Type'] != 'list':
            new_meta = pd.concat([new_meta, pd.DataFrame([meta.iloc[i]])])
        for k in range(duplicate_count):
            if template[i] == list_params[k]:
                indices = get_duplicate_indices(template, list_params[k])
                for l in range(list_counts[k]):
                    sub_index = (f'[{l}]')
                    for m in indices:
                        new_row = meta.iloc[m].copy(deep=True)
                        new_row['Value'] = \
                            new_row['Value'].replace('[i]',sub_index)
                        new_meta = pd.concat([new_meta,
                            pd.DataFrame([new_row])])
    new_meta.reset_index(drop=True, inplace=True)

    return new_meta


def iso_template2lists(excel_file, worksheet):
    """Convert ISO template to lists"""

    ### ISO metadata structure
    meta = pd.read_excel(excel_file, sheet_name=worksheet)
    (row_count, level_count) = meta.shape
    meta_type = list(meta['Type'])
    del meta['Type']
    meta_value = meta['Value']
    del meta['Value']
    if worksheet == 'ISO Metadata Structure':
        iso_attributes = list(meta['Attribute'])
        del meta['Attribute']
        iso_attribute_values = list(meta['AttributeValue'])
        del meta['AttributeValue']
    else:
        iso_attributes = []
        for i in range(row_count):
            iso_attributes.append('nan')
    (row_count, level_count) = meta.shape

    ### Read all entries in the template
    iso_template = []
    for i in range(row_count):
        param = ''
        for k in range(level_count):
            element = str(meta.at[i, k])
            if element != 'nan':
                param += element
                if k < (level_count-1) and str(meta.at[i, k+1]) != 'nan':
                    param += '/'
        iso_template.append(param)

    ### Collect elements that contain lists
    duplicates = []
    for i in range(level_count):
        for k in range(row_count):
            occurences = [x.start() for x in re.finditer('/', iso_template[k])]
            if len(occurences) == i:
                if 'list' in meta_type[k]:
                    duplicates.append(iso_template[k])
        list_params = [item for item, \
            count in collections.Counter(duplicates).items() if count > 0]
        if len(list_params) > 0:
            iso_template = update_template(iso_template, list_params)

    ### Read parameters out of template
    iso_params = []
    iso_values = []
    for i in range(row_count):
        if 'attribute' in meta_type[i]:
            attrib = iso_template[i] + '/' + str(iso_attributes[i])
            iso_params.append(attrib)
            iso_values.append(iso_attribute_values[i])
        if 'value' in meta_type[i]:
            iso_params.append(iso_template[i])
            iso_values.append(meta_value[i])

    return (iso_template, iso_params, iso_values)


def merge_iso_dem_template(iso_template, dem_template):
    """Merge ISO and DEM templates"""

    ### Analyze ISO template
    max_identification_info = 0
    max_identification_index = 0
    max_content_info = 0
    max_content_index = 0
    iso_template_count = len(iso_template)
    for i in range(iso_template_count):
        test_level = iso_template[i][-2:-1]
        if iso_template[i][:-3].endswith('identificationInfo'):
            if max_identification_info < int(test_level):
                max_identification_info = int(test_level)
        if 'identificationInfo' in iso_template[i] and \
            max_identification_index < i:
            max_identification_index = i
        if iso_template[i][:-3].endswith('contentInfo'):
            if max_content_info < int(test_level):
                max_content_info = int(test_level)
        if 'contentInfo' in iso_template[i] and max_content_index < i:
            max_content_index = i

    ### Analyze DEM template
    max_index = 0
    dem_template_count = len(dem_template)
    for i in range(dem_template_count):
        if 'identificationInfo' in dem_template[i] and max_index < i:
            max_index = i

    ### Combine ISO and DEM templates
    merged_template = []
    for i in range(max_identification_index+1):
        merged_template.append(iso_template[i])
    original_string = 'identificationInfo[1]'
    replace_string = (f'identificationInfo[{max_identification_info+1}]')
    for i in range(max_index+1):
        if 'identificationInfo' in dem_template[i]:
            dem_info = dem_template[i].replace(original_string, replace_string)
            merged_template.append(dem_info)
    for i in range(max_identification_index+1, max_content_index+1):
        merged_template.append(iso_template[i])
    original_string = 'contentInfo[1]'
    replace_string = (f'contentInfo[{max_identification_info+1}]')
    for i in range(max_index+1, len(dem_template)):
        dem_info = dem_template[i].replace(original_string, replace_string)
        merged_template.append(dem_info)
    for i in range(max_content_index+1, len(iso_template)):
        merged_template.append(iso_template[i])

    return merged_template


def add_dem_lists(iso_template, iso_params, iso_values, dem_template,
    dem_params, dem_values):
    """Add DEM lists"""

    ### Analyze ISO template
    max_identification_info = 0
    max_identification_index = 0
    max_content_info = 0
    max_content_index = 0
    iso_template_count = len(iso_template)
    for i in range(iso_template_count):
        test_level = iso_template[i][-2:-1]
        if iso_template[i][:-3].endswith('identificationInfo'):
            if max_identification_info < int(test_level):
                max_identification_info = int(test_level)
        if 'identificationInfo' in iso_template[i] and \
            max_identification_index < i:
            max_identification_index = i
        if iso_template[i][:-3].endswith('contentInfo'):
            if max_content_info < int(test_level):
                max_content_info = int(test_level)
        if 'contentInfo' in iso_template[i] and max_content_index < i:
            max_content_index = i

    ### Analyze DEM template
    max_index = 0
    dem_template_count = len(dem_template)
    for i in range(dem_template_count):
        if 'identificationInfo' in dem_template[i] and max_index < i:
            max_index = i

    ### Combine ISO and DEM templates
    merged_template = []
    for i in range(max_identification_index+1):
        merged_template.append(iso_template[i])
    original_string = 'identificationInfo[1]'
    replace_string = (f'identificationInfo[{max_identification_info+1}]')
    for i in range(max_index+1):
        if 'identificationInfo' in dem_template[i]:
            dem_info = dem_template[i].replace(original_string, replace_string)
            merged_template.append(dem_info)
    for i in range(max_identification_index+1, max_content_index+1):
        merged_template.append(iso_template[i])
    original_string = 'contentInfo[1]'
    replace_string = (f'contentInfo[{max_identification_info+1}]')
    for i in range(max_index+1, len(dem_template)):
        dem_info = dem_template[i].replace(original_string, replace_string)
        merged_template.append(dem_info)
    for i in range(max_content_index+1, len(iso_template)):
        merged_template.append(iso_template[i])

    ### Analyze ISO params and values
    max_identification_index = 0
    max_content_index = 0
    iso_param_count = len(iso_params)
    for i in range(iso_param_count):
        if 'identificationInfo' in iso_params[i] and \
            max_identification_index < i:
            max_identification_index = i
        if 'contentInfo' in iso_params[i] and max_content_index < i:
            max_content_index = i

    ### Analyze DEM params
    max_index = 0
    dem_param_count = len(dem_params)
    for i in range(dem_param_count):
        if 'identificationInfo' in dem_params[i] and max_index < i:
            max_index = i

    ### Combine ISO and DEM params and values
    merged_params = []
    merged_values = []
    for i in range(max_identification_index+1):
        merged_params.append(iso_params[i])
        merged_values.append(iso_values[i])
    original_string = 'identificationInfo[1]'
    replace_string = (f'identificationInfo[{max_identification_info+1}]')
    for i in range(max_index+1):
        dem_info = dem_params[i].replace(original_string, replace_string)
        merged_params.append(dem_info)
        merged_values.append(dem_values[i])
    for i in range(max_identification_index+1, max_content_index+1):
        merged_params.append(iso_params[i])
        merged_values.append(iso_values[i])
    original_string = 'contentInfo[1]'
    replace_string = (f'contentInfo[{max_identification_info+1}]')
    for i in range(max_index+1, len(dem_params)):
        dem_info = dem_params[i].replace(original_string, replace_string)
        merged_params.append(dem_info)
        merged_values.append(dem_values[i])
    for i in range(max_content_index+1, len(iso_params)):
        merged_params.append(iso_params[i])
        merged_values.append(iso_values[i])

    return (merged_template, merged_params, merged_values)


def iso_attributes2lists(excel_file):
    """Convert ISO attributes to lists"""

    ### ISO metadata structure
    meta = pd.read_excel(excel_file, sheet_name='ISO Metadata Structure')
    iso_attributes = list(meta['Attribute'])
    iso_attribute_values = list(meta['AttributeValue'])
    (iso_template, _, _) = \
        iso_template2lists(excel_file, 'ISO Metadata Structure')

    return (iso_attributes, iso_attribute_values, iso_template)


def iso_xml_structure(excel_file, iso_template, iso_params, iso_values,
    ns_flag):
    """Build ISO XML structure"""

    ### Get name_space look up table
    ns_lut = pd.read_excel(excel_file, sheet_name='Namespaces LUT')
    ns_pairs = dict(zip(ns_lut['Parameter'], ns_lut['Namespace']))

    ### Build name_space
    ns_lut = pd.read_excel(excel_file, sheet_name='Namespaces URL')
    ns_lut = ns_lut.set_index('Namespace')
    ns_xsi = {'xsi': ns_lut.at['xsi','URL']}
    ns_gmd = {'gmd': ns_lut.at['gmd','URL']}
    ns_gco = {'gco': ns_lut.at['gco','URL']}
    ns_xs = {'xs': ns_lut.at['xs','URL']}
    ns_eos = {'eos': ns_lut.at['eos','URL']}
    ns_echo = {'echo': ns_lut.at['echo','URL']}
    ns_xlink = {'xlink': ns_lut.at['xlink','URL']}
    ns_gml = {'gml': ns_lut.at['gml','URL']}
    ns_gmi = {'gmi': ns_lut.at['gmi','URL']}
    ns_gmx = {'gmx': ns_lut.at['gmx','URL']}
    ns_all = dict(
        list(ns_xsi.items()) +
        list(ns_gmd.items()) +
        list(ns_gco.items()) +
        list(ns_xs.items()) +
        list(ns_eos.items()) +
        list(ns_echo.items()) +
        list(ns_xlink.items()) +
        list(ns_gml.items()) +
        list(ns_gmi.items()) +
        list(ns_gmx.items())
    )

    iso_elements = []
    iso_element_count = []
    row_count = len(iso_template)
    level_count = 0
    for i in range(row_count):
        meta_elements = iso_template[i].split('/')
        if len(meta_elements) > level_count:
            level_count = len(meta_elements)
        iso_elements.append(meta_elements)
        iso_element_count.append(len(meta_elements))
    level_elements = []
    for k in range(level_count):
        meta = []
        for i in range(row_count):
            if k < len(iso_elements[i]):
                meta.append(iso_elements[i][k])
        level_elements.append(unique_list(meta))

    ### Build XMl element tree from nested dictionary
    name_space = (f"{ns_lut.at[ns_pairs['DS_Series'],'URL']}")
    if ns_flag:
        root = \
            et.Element(f"{name_space}{'DS_Series'}", nsmap=ns_all)
    else:
        root = et.Element('DS_Series')

    ## Build structure
    ref = {}
    ref['DS_Series'] = root
    for i in range(1, level_count):
        indices = [k for k, x in enumerate(iso_element_count) if x == i+1]
        index_count = len(indices)
        for k in range(index_count):
            element = iso_template[indices[k]]
            parent = ref[element.rsplit('/',1)[0]]
            param = element.rsplit('/',1)[1]
            if 'EOS_AdditionalAttributes' in element:
                name_space = (f"{ns_lut.at['eos','URL']}")
            elif ('acquisitionInformation' in element) and \
                param in ('identifier','description','type'):
                name_space = (f"{ns_lut.at['gmi','URL']}")
            else:
                name_space = \
                    (f"{ns_lut.at[ns_pairs[param.split('[')[0]],'URL']}")
            if ns_flag:
                child = et.SubElement(parent,
                    f"{name_space}{param.split('[')[0]}", ns_map=ns_all)
            else:
                child = et.SubElement(parent, param.split('[')[0])
            ref[element] = child

    ## Add values
    iso_param_count = len(iso_params)
    for i in range(iso_param_count):
        element = iso_params[i]
        param = element.rsplit('/',1)[1]
        if (param == 'value') and (('result/value' not in element) or \
            ('EOS_AdditionalAttribute/value' not in element)):
            continue
        elif param == 'id':
            reference = ref[element.rsplit('/',1)[0]]
            name_space = \
                (f"{ns_lut.at[ns_pairs[param.split('[')[0]],'URL']}")
            if ns_flag:
                ns_param = (f"{name_space}{param.split('[')[0]}",
                    "ns_map=ns_all")
            else:
                ns_param = param.split('[')[0]
            reference.attrib[ns_param] = str(iso_values[i])
        elif param in ('href','nilReason'):
            reference = ref[element.rsplit('/',1)[0]]
            name_space = \
                (f"{ns_lut.at[ns_pairs[param.split('[')[0]],'URL']}")
            if ns_flag:
                ns_param = (f"{name_space}{param.split('[')[0]}",
                    "ns_map=ns_all")
            else:
                ns_param = param.split('[')[0]
            reference.attrib[ns_param] = str(iso_values[i])
        elif param == 'type':
            reference = ref[element.rsplit('/',1)[0]]
            name_space = (f"{ns_lut.at['xsi','URL']}")
            if ns_flag:
                ns_param = (f"{name_space}{param.split('[')[0]}",
                    "ns_map=ns_all")
            else:
                ns_param = param.split('[')[0]
            reference.attrib[ns_param] = str(iso_values[i])
        elif param == 'uom':
            reference = ref[element.rsplit('/',1)[0]]
            reference.attrib[param] = str(iso_values[i])
            reference.text = str(iso_values[i-1])
        elif param.split('[')[0] in ns_pairs:
            name_space = \
                (f"{ns_lut.at[ns_pairs[param.split('[')[0]],'URL']}")
            reference = ref[element]
            reference.text = str(iso_values[i])
            if 'Code' in param:
                if 'EOS' in element:
                    code_url = 'http://earthdata.nasa.gov/metadata/' \
                        'resources/codelist.xml#' + param
                else:
                    code_url = 'https://cdn.earthdata.nasa.gov/iso/' \
                        'resources/Codelist/gmxCodelists.xml#' + param
                reference.attrib['codeList'] = code_url
                reference.attrib['codeListValue'] = str(iso_values[i])
            if param == 'RecordType':
                name_space = (f"{ns_lut.at['xlink','URL']}")
                if ns_flag:
                    ns_param = (f"{name_space}{'href'}",
                        "ns_map=ns_all")
                else:
                    ns_param = param.split('[')[0]
                reference.attrib[ns_param] = "http://earthdata.nasa.gov/" \
                    "schemas/eos/eos.xsd#xpointer(//" \
                    "element[@name='EOS_AdditionalAttributes'])"

    return root


def add_dem_attributes(iso_attributes, dem_attributes, iso_template,
    dem_template):
    """Add DEM attributes"""

    ### Analyze ISO template
    max_identification_info = 0
    max_identification_index = 0
    max_content_info = 0
    max_content_index = 0
    iso_template_count = len(iso_template)
    for i in range(iso_template_count):
        test_level = iso_template[i][-2:-1]
        if iso_template[i][:-3].endswith('identificationInfo'):
            if max_identification_info < int(test_level):
                max_identification_info = int(test_level)
        if 'identificationInfo' in iso_template[i] and \
            max_identification_index < i:
            max_identification_index = i
        if iso_template[i][:-3].endswith('contentInfo'):
            if max_content_info < int(test_level):
                max_content_info = int(test_level)
        if 'contentInfo' in iso_template[i] and max_content_index < i:
            max_content_index = i

    ### Analyze DEM template
    max_index = 0
    dem_template_count = len(dem_template)
    for i in range(dem_template_count):
        if 'identificationInfo' in dem_template[i] and max_index < i:
            max_index = i

    ### Combine ISO and DEM attributes
    merged_attributes = []
    for i in range(max_identification_index+1):
        merged_attributes.append(iso_attributes[i])
    for i in range(max_index+1):
        merged_attributes.append(dem_attributes[i])
    for i in range(max_identification_index+1, max_content_index+1):
        merged_attributes.append(iso_attributes[i])
    for i in range(max_index+1, len(dem_attributes)):
        merged_attributes.append(dem_attributes[i])
    for i in range(max_content_index+1, len(iso_attributes)):
        merged_attributes.append(iso_attributes[i])

    return merged_attributes


def iso_dictionary_structure(excel_file, dem_file, meta_template, meta_params,
  meta_values):
    """Build ISO dictionary structure"""

    ### Read attributes and their values from Excel spreadsheet
    (iso_attributes, _, iso_template) = iso_attributes2lists(excel_file)
    if dem_file:
        (dem_attributes, _, dem_template) = iso_attributes2lists(dem_file)
        iso_attributes = add_dem_attributes(iso_attributes, dem_attributes[5:],
            iso_template, dem_template[5:])

    ### Build dataframe with parameters
    df_params = get_params_dataframe(meta_params)
    (_, old_level) = df_params.shape
    old_level -= 1
    dict_list = \
        params2dict_list(df_params, meta_params, meta_values, old_level)

    ### Build dictionary structure - one level at a time
    for new_level in reversed(range(old_level)):

        new_dict_list = \
            params2dict_list(df_params, meta_params, meta_values, new_level)
        dict_list += new_dict_list
        dict_list = upgrade_dictionary2level(dict_list, old_level, new_level)
        level_params = \
            get_level_params_list(meta_template, iso_attributes, new_level)
        dict_list = merge_dictionary_params(dict_list, level_params, old_level,
            new_level)
        old_level = new_level

    iso_dict = {}
    iso_dict[dict_list[0]['key']] = dict_list[0]['dict']

    return iso_dict


def umm_dictionary_structure(meta_template, meta_params, meta_values):
    """Build UMM dictionary structure"""

    ### Build dataframe with parameters
    df_params = get_params_dataframe(meta_params)
    (_, old_level) = df_params.shape
    old_level -= 1
    dict_list = \
        params2dict_list(df_params, meta_params, meta_values, old_level)

    ### Build dictionary structure - one level at a time
    for new_level in reversed(range(old_level)):

        new_dict_list = \
            params2dict_list(df_params, meta_params, meta_values, new_level)
        dict_list += new_dict_list
        dict_list = upgrade_dictionary2level(dict_list, old_level, new_level)
        level_params = get_level_params_list(meta_template, None, new_level)
        dict_list = merge_dictionary_params(dict_list, level_params, old_level,
            new_level)
        old_level = new_level

    iso_dict = {}
    iso_dict[dict_list[0]['key']] = dict_list[0]['dict']

    return iso_dict


def meta_json_file(meta_structure, json_file):
    """Write metadata structure to JSON file"""

    ### Write JSON file
    with open(json_file, 'w', encoding='utf-8') as file:
        json.dump(meta_structure, file, indent=2)


def meta_xml_file(meta_structure, xml_file):
    """Write metaadata structure to XML File"""

    ### Write XML file
    with open(xml_file, 'wb', encoding='utf-8') as file:
        file.write(et.tostring(meta_structure,
            xml_declaration=True, encoding='utf-8',
            pretty_print=True))


def generate_product_dictionary(product_type, data_type, meta_file,
    log_file):
    """Generate product dictionary"""

    if product_type == 'gamma_rtc':
        meta = gamma_rtc_log2meta(data_type, meta_file, log_file)

    return meta


def get_metadata_values(meta_file, params):
    """Get metadata values"""

    product_params = []

    parser = et.XMLParser(remove_blank_text=True)
    tree = et.parse(meta_file, parser)
    root = tree.getroot()
    for elem in root.getiterator():
        if len(elem.attrib) > 0:
            key = elem.attrib.keys()
            item = elem.attrib.items()
            (key, value) = item[0]
            if '}' in key:
                del elem.attrib[key]
                key = key.split('}')[1]
                elem.set(key, value)
        if not hasattr(elem.tag, 'find'):
            continue
        i = elem.tag.find('}')
        if i >= 0:
            elem.tag = elem.tag[i+1:]
    objectify.deannotate(root, cleanup_name_spaces=True)

    param_count = len(params)
    for i in range(param_count):
        try:
            param = tree.xpath('/'+params[i])[0].text
            product_params.append(get_value(param))
        except:
            (element, attribute) = params[i].rsplit('/',1)
            if attribute != 'value':
                element_attribute = ('/'+element+'/@'+attribute)
                param = tree.xpath(element_attribute)[0]
                product_params.append(get_value(param))
            else:
                param = tree.xpath('/'+element)[0].text
                product_params.append(get_value(param))

    return product_params


def product_dictionary2values(values, prod_dict):
    """Convert product dictionary to values"""

    dict_keys = list(prod_dict.keys())
    prod_values = []
    value_count = len(values)
    for i in range(value_count):
        keys = []
        dict_key_count = len(dict_keys)
        for k in range(dict_key_count):
            if str(values[i]).count(dict_keys[k]) > 0:
                keys.append(dict_keys[k])
        for key in keys:
            start = values[i].index('{# ')
            stop = values[i].index(' #}')
            parameter = values[i][start:stop+3]
            values[i] = \
                values[i].replace(parameter, str(prod_dict[key]))
            try:
                if isinstance(prod_dict[key], float):
                    values[i] = float(values[i])
                elif isinstance(prod_dict[key], bool):
                    values[i] = bool(values[i])
                elif isinstance(prod_dict[key], int):
                    values[i] = int(values[i])
            except ValueError:
                pass

        prod_values.append(values[i])

    return prod_values


def granule_iso2umm_values(iso_params, iso_values, iso2umm, umm_values):
    """Convert granule ISO to UMM values"""

    (param_count, _) = iso2umm.shape
    umm_value_count = len(umm_values)
    for k in range(umm_value_count):
        for i in range(param_count):
            umm_param_string = (f"{{# UMM_{iso2umm.at[i, 'UMM']} #}}")
            if umm_param_string == umm_values[k]:
                iso_index = iso_params.index(iso2umm.at[i, 'ISO'])
                umm_values[k] = iso_values[iso_index]

    return umm_values


def properties2values(template, properties):
    """Convert properties to values"""

    values = []
    for param in template:
        if param in properties.keys():
            try:
                if isinstance(properties[param], float):
                    values.append(float(properties[param]))
                elif isinstance(properties[param], bool):
                    values.append(bool(properties[param]))
                elif isinstance(properties[param], int):
                    values.append(int(properties[param]))
                else:
                    values.append(properties[param])
            except ValueError:
                pass
        else:
            values.append(param)

    return values


def clean_json_structure(iso_structure):
    """Clean JSON structure"""

    meta = iso_structure['DS_Series']['composedOf']['DS_DataSet']['has'] \
        ['MI_Metadata']

    identification_info = meta['identificationInfo']
    new_identification = []
    identification_info_count = len(identification_info)
    for i in range(identification_info_count):
        value = identification_info[i]['MD_DataIdentification']['citation'] \
            ['CI_Citation']['title']['FileName']
        if not value.startswith('{# ISO'):
            new_identification.append(identification_info[i])
    meta['identificationInfo'] = new_identification

    content_info = meta['contentInfo']
    new_content = []
    content_info_count = len(content_info)
    for i in range(content_info_count):
        mi_band = \
            content_info[i]['MD_CoverageDescription']['dimension']['MI_Band']
        if 'maxValue' in mi_band:
            value = mi_band['maxValue']['Real']
        elif 'otherPropertyType' in mi_band:
            value = mi_band['otherProperty']['Record']['additionalAttributes'] \
                ['EOS_AdditionalAttributes']['additionalAttribute'][0] \
                ['EOS_AdditionalAttribute']['value']['CharacterString']
        if isinstance(value, float):
            new_content.append(content_info[i])
    meta['contentInfo'] = new_content

    return iso_structure


def clean_xml_structure(iso_file):
    """Clean XML structure"""

    ns_gmd = {'gmd' : 'http://www.isotc211.org/2005/gmd'}
    ns_gco = {'gco' : 'http://www.isotc211.org/2005/gco'}
    ns_eos = {'eos' : 'http://earthdata.nasa.gov/schema/eos'}
    ns_gmi = {'gmi' : 'http://www.isotc211.org/2005/gmi'}
    ns_gmx = {'gmx' : 'http://www.isotc211.org/2005/gmx'}
    ns_all = dict(
        list(ns_gmd.items()) +
        list(ns_gco.items()) +
        list(ns_eos.items()) +
        list(ns_gmi.items()) +
        list(ns_gmx.items())
    )

    parser = et.XMLParser(remove_blank_text=True)
    doc = et.parse(iso_file, parser)
    meta = doc.xpath('/gmd:DS_Series/gmd:composedOf/gmd:DS_DataSet' \
        '/gmd:has/gmi:MI_Metadata', namespaces=ns_all)[0]
    identification_info_count = 0
    content_info_count = 0
    remove = []
    meta_count = len(meta)
    for i in range(meta_count):
        if 'identificationInfo' in meta[i].tag:
            xml = ('/gmd:DS_Series/gmd:composedOf/'
                'gmd:DS_DataSet/gmd:has/gmi:MI_Metadata/'
                f'gmd:identificationInfo[{identification_info_count+1}]/'
                'gmd:MD_DataIdentification/'
                'gmd:citation/gmd:CI_Citation/gmd:title/gmx:FileName')
            file_name = doc.xpath(xml, namespaces=ns_all)[0].text
            if file_name.startswith('{# ISO'):
                remove.append(i)
            identification_info_count += 1
        elif 'contentInfo' in meta[i].tag:
            xml = ('/gmd:DS_Series/gmd:composedOf/'
                'gmd:DS_DataSet/gmd:has/gmi:MI_Metadata/'
                f'gmd:contentInfo[{content_info_count+1}]/'
                'gmd:MD_CoverageDescription/gmd:dimension/gmi:MI_Band')
            mi_band = doc.xpath(xml, namespaces=ns_all)[0]
            if 'maxValue' in mi_band[0].tag:
                xml = ('/gmd:DS_Series/gmd:composedOf/'
                    'gmd:DS_DataSet/gmd:has/gmi:MI_Metadata/'
                    f'gmd:contentInfo[{content_info_count+1}]/'
                    'gmd:MD_CoverageDescription/gmd:dimension/gmi:MI_Band/'
                    'gmd:maxValue/gco:Real')
                value = doc.xpath(xml, namespaces=ns_all)[0].text
            elif 'otherPropertyType' in mi_band[0].tag:
                xml = ('/gmd:DS_Series/gmd:composedOf/'
                    'gmd:DS_DataSet/gmd:has/gmi:MI_Metadata/'
                    f'gmd:contentInfo[{content_info_count+1}]/'
                    'gmd:MD_CoverageDescription/gmd:dimension/gmi:MI_Band/'
                    'eos:otherProperty/gco:Record/eos:additionalAttributes/'
                    'eos:EOS_AdditionalAttributes/eos:additionalAttribute[1]/'
                    'eos:EOS_AdditionalAttribute/eos:value/'
                    'eos:CharacterString')
                value = doc.xpath(xml, namespaces=ns_all)[0].text
            try:
                float(value)
            except ValueError:
                remove.append(i)
            content_info_count += 1

    remove_count = len(remove)-1
    for i in range(remove_count,-1,-1):
        del meta[remove[i]]

    return doc


def gamma_rtc_log2meta(data_type, meta_file, log_file):
    """Convert GAMMA RTC log to metadata structure"""

    ns_safe = {'safe': 'http://www.esa.int/safe/sentinel-1.0'}
    ns_all = dict(
        list(ns_safe.items())
    )

    ### Read information from log file
    with open(log_file, encoding='utf-8') as file:
        lines = file.readlines()

    meta = {}
    meta['ISO_RTC_metadataCreationTime'] = datetime.utcnow().isoformat() + 'Z'
    meta['ISO_RTC_productCreationTime'] = parse_line(lines[0], 'dateString')
    meta['ISO_RTC_radiometry'] = parse_line(lines[3], 'value')
    meta['ISO_RTC_scale'] = parse_line(lines[4], 'value')
    meta['ISO_RTC_DEM_type'] = parse_line(lines[13], 'value').upper()
    for line in lines:
        value = parse_line(line, 'value')
        if value and value.startswith('geoid file'):
            meta['ISO_RTC_GammaVersion'] = \
                [i for i in value.split('/') \
                    if 'GAMMA_SOFTWARE' in i][0].split('-')[1]
        elif value and value.startswith('number of azimuth looks'):
            meta['ISO_RTC_azimuthLooks'] = int(value.split(':')[1])
        elif value and value.startswith('number of range looks'):
            meta['ISO_RTC_rangeLooks'] = int(value.split(':')[1])
    meta['ISO_RTC_XML_filename'] = log_file.replace('.log','.iso.xml')

    if data_type == 'sentinel':

        ### Extract more metadata from manifest file
        input_path = os.path.dirname(meta_file)
        meta['ISO_RTC_inputGranule'] = os.path.basename(input_path)[:-5]
        manifest_file = meta_file
        manifest = read_manifest_file(manifest_file)
        doc = et.fromstring(manifest['metadataSection'])
        meta['ISO_RTC_startTime'] = \
            doc.xpath('/metadataSection/metadataObject' \
            '[@ID="acquisitionPeriod"]/metadataWrap/xmlData/' \
            'safe:acquisitionPeriod/safe:startTime',
            namespaces=ns_all)[0].text + 'Z'
        meta['ISO_RTC_stopTime'] = \
            doc.xpath('/metadataSection/metadataObject' \
            '[@ID="acquisitionPeriod"]/metadataWrap/xmlData/' \
            'safe:acquisitionPeriod/safe:stopTime',
            namespaces=ns_all)[0].text + 'Z'

        ### Extract more metadata from annotation file
        annotation_path = os.path.join(input_path, 'annotation', '*.xml')
        annotation_file = glob.glob(f'{annotation_path}')[0]
        annotation = read_annotation_file(annotation_file)
        doc = et.fromstring(annotation['adsHeader'])
        meta['ISO_RTC_mission'] = \
            doc.xpath('/adsHeader/missionId')[0].text.replace('S','SENTINEL-')
        meta['ISO_RTC_beamMode'] = doc.xpath('/adsHeader/mode')[0].text
        meta['ISO_RTC_absoluteOrbitNumber'] = \
            int(doc.xpath('/adsHeader/absoluteOrbitNumber')[0].text)
        doc = et.fromstring(annotation['generalAnnotation'])
        meta['ISO_RTC_orbitPassDirection'] = \
            doc.xpath('/generalAnnotation/productInformation/pass')[0].text.lower()
        speed_of_light = 299792458.0
        frequency = float(doc.xpath('/generalAnnotation/productInformation/' \
            'radarFrequency')[0].text)
        meta['ISO_RTC_wavelength'] = speed_of_light / frequency

        ### Examine product directory for various files and calculate stats if exists
        polarization = meta['ISO_RTC_inputGranule'][14:16]
        if polarization in ('SH','DH'):
            main_pol = 'HH'
            pol_string = 'horizontal'
        elif polarization in ('SV','DV'):
            main_pol = 'VV'
            pol_string = 'vertical'

    ### Extract Hyp3 version out of README file
    readme_file = log_file.replace('.log','.README.md.txt')
    with open(readme_file, encoding='utf-8') as file:
        lines = file.readlines()
    for line in lines:
        if parse_line(line, 'key') == 'Metadata version':
            meta['ISO_RTC_productVersion'] = parse_line(line, 'value')

    ### Setting additional metadata
    meta['ISO_RTC_lookDirection'] = 'RIGHT'
    if meta['ISO_RTC_mission'] == 'ALOS':
        meta['ISO_RTC_sensor'] = 'PALSAR'
    else:
        meta['ISO_RTC_sensor'] = 'SAR'

    ## Look at product file
    product_file = log_file.replace('.log',f'_{main_pol}.tif')
    if os.path.isfile(product_file):
        meta['ISO_RTC_terrainCorrectedImage'] = os.path.basename(product_file)
        raster = gdal.Open(product_file)
        meta['ISO_RTC_azimuthCount'] = raster.RasterXSize
        band = raster.GetRasterBand(1)
        stats = band.GetStatistics(False, True)
        proj = osr.SpatialReference()
        proj.ImportFromWkt(raster.GetProjectionRef())
        if proj.GetAttrValue('AUTHORITY', 0) == 'EPSG':
            epsg = int(proj.GetAttrValue('AUTHORITY', 1))
            geo_transform = raster.GetGeoTransform()
            meta['ISO_RTC_azimuthPixelSize'] = geo_transform[1]
            meta['ISO_RTC_rangePixelSize'] = -geo_transform[5]
            meta['ISO_RTC_azimuthCount'] = raster.RasterXSize
            meta['ISO_RTC_rangeCount'] = raster.RasterYSize
        else:
            epsg = 0
        meta['ISO_RTC_minValue'] = stats[0]
        meta['ISO_RTC_maxValue'] = stats[1]
        meta['ISO_RTC_meanValue'] = stats[2]
        meta['ISO_RTC_standardDeviation'] = stats[3]
        meta['ISO_RTC_transmittedPolarization'] = pol_string
        meta['ISO_RTC_receivedPolarization'] = pol_string
        meta['ISO_RTC_epsgCode'] = epsg

    ## Look at DEM file
    dem_file = log_file.replace('.log','_dem.tif')
    if os.path.isfile(dem_file):
        meta['ISO_RTC_digitalElevationModel'] = os.path.basename(dem_file)
        raster = gdal.Open(product_file)
        band = raster.GetRasterBand(1)
        stats = band.GetStatistics(False, True)
        meta['ISO_RTC_DEM_minValue'] = stats[0]
        meta['ISO_RTC_DEM_maxValue'] = stats[1]
        meta['ISO_RTC_DEM_meanValue'] = stats[2]
        meta['ISO_RTC_DEM_standardDeviation'] = stats[3]

    ## Look at the incidence angle map
    incidence_angle_file = log_file.replace('.log','_inc.tif')
    if os.path.isfile(incidence_angle_file):
        meta['ISO_RTC_incidenceAngleMap'] = \
            os.path.basename(incidence_angle_file)
        raster = gdal.Open(product_file)
        band = raster.GetRasterBand(1)
        stats = band.GetStatistics(False, True)
        meta['ISO_RTC_incidenceMinValue'] = stats[0]
        meta['ISO_RTC_incidenceMaxValue'] = stats[1]
        meta['ISO_RTC_incidenceMeanValue'] = stats[2]
        meta['ISO_RTC_incidenceStandardDeviation'] = stats[3]

    ## Look at the layover shadow mask
    layover_shadow_mask_file = log_file.replace('.log','_ls_map.tif')
    if os.path.isfile(layover_shadow_mask_file):
        meta['ISO_RTC_layoverShadowMask'] = \
            os.path.basename(layover_shadow_mask_file)
        raster = gdal.Open(layover_shadow_mask_file)
        band = raster.GetRasterBand(1).ReadAsArray()
        pixel_count = raster.RasterXSize * raster.RasterYSize
        (histogram, _) = np.histogram(band, bins=np.arange(256))
        no_layover_shadow = 0.0
        true_layover = 0.0
        layover = 0.0
        true_shadow = 0.0
        shadow = 0.0
        pixel_count -= float(histogram[0])
        histogram_count = len(histogram)
        for i in range(histogram_count):
            if i & 1:
                no_layover_shadow += float(histogram[i]) / pixel_count
            if i & 2:
                true_layover += float(histogram[i]) / pixel_count
            if i & 4:
                layover += float(histogram[i]) / pixel_count
            if i & 8:
                true_shadow += float(histogram[i]) / pixel_count
            if i & 16:
                shadow += float(histogram[i]) / pixel_count
        meta['ISO_RTC_noLayoverShadowPercentage'] = no_layover_shadow * 100.0
        meta['ISO_RTC_trueLayoverPercentage'] = true_layover * 100.0
        meta['ISO_RTC_layoverPercentage'] = layover * 100.0
        meta['ISO_RTC_trueShadowPercentage'] = true_shadow * 100.0
        meta['ISO_RTC_shadowPercentage'] = shadow * 100.0

    ## Look at the scattering area map
    scattering_area_file = log_file.replace('.log','_area.tif')
    if os.path.isfile(scattering_area_file):
        meta['ISO_RTC_scatteringAreaMap'] = \
            os.path.basename(scattering_area_file)
        raster = gdal.Open(product_file)
        band = raster.GetRasterBand(1)
        stats = band.GetStatistics(False, True)
        meta['ISO_RTC_scatteringMinValue'] = stats[0]
        meta['ISO_RTC_scatteringMaxValue'] = stats[1]
        meta['ISO_RTC_scatteringMeanValue'] = stats[2]
        meta['ISO_RTC_scatteringStandardDeviation'] = stats[3]

    ## Browse
    if polarization in ('DH','DV'):
        meta['ISO_RTC_browseImage'] = \
            os.path.basename(log_file.replace('.log','_rgb.png'))
        meta['ISO_RTC_KMLoverlay'] = \
            os.path.basename(log_file.replace('.log','_rgb.kmz'))
    else:
        meta['ISO_RTC_browseImage'] = \
            os.path.basename(log_file.replace('.log','.png'))
        meta['ISO_RTC_KMLoverlay'] = \
            os.path.basename(log_file.replace('.log','.kmz'))

    ## Extent
    (lon_min, lon_max, lat_min, lat_max) = get_latlon_extent(product_file)
    meta['ISO_RTC_westBoundLongitude'] = lon_min
    meta['ISO_RTC_eastBoundLongitude'] = lon_max
    meta['ISO_RTC_northBoundLatitude'] = lat_max
    meta['ISO_RTC_southBoundLatitude'] = lat_min

    return meta


def project_template2dict(excel_file, worksheet, level, data):
    """Convert project template to JSON dictionary"""

    ### Project metadata structure
    print(f'Metadata level: {level}')
    meta = pd.read_excel(excel_file, sheet_name=worksheet)
    if level == 'short':
        meta = meta[meta['Level']=='short']
        meta.reset_index(inplace=True, drop=True)

    (row_count, column_count) = meta.shape
    level_count = column_count - 6

    ### Read all entries from the dataframe into the template
    project_template = dataframe2template(meta, level_count)

    ### Collect elements that contain lists
    duplicates = []
    list_counts = []
    flat_data = flatdict.FlatDict(data, delimiter='.')
    for i in range(level_count):
        for k in range(row_count):
            occurences = \
                [x.start() for x in re.finditer('/', project_template[k])]
            if len(occurences) == i:
                if 'list' in meta.at[k,'Type']:
                    list_counts.append(int(flat_data[meta.at[k,'Value'][3:-3]]))
                    duplicates.append(project_template[k])
    list_params = [item for item, \
        count in collections.Counter(duplicates).items() if count > 0]
    if len(list_params) > 0:
        meta = add_duplicates2metadata(meta, list_params, list_counts)
        (row_count, _) = meta.shape
        project_template = dataframe2template(meta, level_count)
        project_template = update_template(project_template, list_params)

    ### Read parameters out of template
    project_params = []
    project_values = []
    for i in range(row_count):
        if 'attribute' in meta.at[i,'Type']:
            attrib = project_template[i] + '/@' + str(meta.at[i,'Attribute'])
            project_params.append(attrib)
            project_values.append(meta.at[i,'AttributeValue'])
        if 'value' in meta.at[i,'Type']:
            project_params.append(project_template[i])
            project_values.append(meta.at[i,'Value'])

    ### Build dataframe with parameters
    df_params = get_params_dataframe(project_params)
    (_, old_level) = df_params.shape
    old_level -= 1
    dict_list = \
        params2dict_list(df_params, project_params, project_values, old_level)

    ### Build dictionary structure - one level at a time
    for new_level in reversed(range(old_level)):

        new_dict_list = params2dict_list(df_params, project_params,
            project_values, new_level)
        dict_list += new_dict_list
        dict_list = upgrade_dictionary2level(dict_list, old_level, new_level)
        level_params = get_level_params_list(project_template,
            list(meta['Attribute']), new_level)
        dict_list = merge_dictionary_params(dict_list, level_params, old_level,
            new_level)
        old_level = new_level

    project_dict = {}
    project_dict[dict_list[0]['key']] = dict_list[0]['dict']

    return project_dict
