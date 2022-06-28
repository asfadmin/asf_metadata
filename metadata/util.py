#!/usr/bin/python3

"""Import packages"""
from itertools import groupby
from operator import itemgetter
import os
import configparser
import pandas as pd
from sqlalchemy import create_engine
import requests
import remotezip


def split_element(key, value, out):
    """Split element"""

    key, *rest = key.split('/', 1)
    if rest:
        split_element(rest[0], value, out.setdefault(key, {}))
    else:
        out[key] = value


def nest_dict(flat_dict):
    """Build nested dictionary"""

    result = {}
    for key, value in flat_dict.items():
        split_element(key, value, result)

    return result


def is_float(value):
    """Is value float?"""

    try:
        float(value)
        return True
    except ValueError:
        return False


def is_int(value):
    """Is value integer?"""

    try:
        int(value)
        return True
    except ValueError:
        return False


def get_value(value_str):
    """Get value"""

    if is_int(value_str):
        out_value = int(value_str)
    elif is_float(value_str):
        out_value = float(value_str)
    else:
        out_value = value_str.strip()
    return out_value


def rreplace(string, find, replace):
    """Replace from right"""

    reversed_string = string[::-1]
    replaced = reversed_string.replace(find[::-1], replace[::-1], 1)
    return replaced[::-1]


def unique_list(str_list):
    """Build unique list"""

    out_list = []
    for element in str_list:
        if element not in out_list:
            out_list.append(element)
    return out_list


def meta2dict(meta_file):
    """Convert metadata to dictionary"""

    ### Extract parameters and values from metadata file
    meta_params = []
    meta_values = []
    i = 0
    section = ''
    with open(meta_file, encoding='utf-8') as file:
        lines = [x.strip() for x in file.readlines()]
    for line in lines:
        if '{' in line:
            group = line.split(' {')[0].lstrip()
            section = group
            if section == 'vector':
                i += 1
        elif ':' in line:
            param = line.split(':')[0].lstrip()
            if section == 'vector':
                meta_params.append(f'state/vector[{i}]/{param}')
            elif section == 'param':
                meta_params.append(f'projection/param/{param}')
            elif section == 'utm':
                meta_params.append(f'projection/param/utm/{param}')
            elif section == 'band_stats':
                meta_params.append(f'statistics/band_stats/{param}')
            else:
                if param == 'meta_version':
                    meta_params.append(param)
                else:
                    meta_params.append(section + '/' + param)
            value = line.split(': ')[1].split('#')[0].rstrip()
            if is_int(value):
                meta_values.append(int(value))
            elif is_float(value):
                meta_values.append(float(value))
            else:
                meta_values.append(value)
        elif '}' in line:
            section = ''
        else:
            continue

    ### Establish flattened dictionary
    flat_dict = dict(zip(meta_params, meta_values))

    ### Convert flattened into nested dictionary
    meta_dict = nest_dict(flat_dict)

    return meta_dict


def merge_dict(first, second):
    """Merge dictionary"""

    key = None
    if first is None or isinstance(first, (str, int, float)):
        first = second
    elif isinstance(first, list):
        if isinstance(second, list):
            if len(first[0]) == len(second[0]):
                first[0] = merge_dict(first[0], second[0])
            elif next(iter(second[0])) not in first[len(first)-1].keys():
                first[len(first)-1] = \
                    merge_dict(first[len(first)-1], second[0])
            else:
                first.extend(second)
        else:
            first.append(second)
    elif isinstance(first, dict):
        if isinstance(second, dict):
            for key in second:
                if key in first:
                    first[key] = merge_dict(first[key], second[key])
                else:
                    first[key] = second[key]

    return first


def get_consecutive_list(meta_list):
    """Get consecutive list"""

    for _, group in groupby(enumerate(meta_list), lambda x: x[0]-x[1]):
        yield list(map(itemgetter(1), group))


def get_level_params_list(meta_template, meta_attributes, level):
    """Get level parameter list"""

    level_count = 0
    row_count = len(meta_template)
    for i in range(row_count):
        meta_element = meta_template[i].split('/')
        level_count = len(meta_element) \
            if len(meta_element) > level_count else level_count

    dataframe = pd.DataFrame(index=range(row_count), columns=range(level_count))
    for i in range(row_count):
        meta_element = meta_template[i].split('/')
        for k in range(len(meta_element)):
            dataframe.at[i,k] = meta_element[k]
        if meta_attributes:
            if isinstance(meta_attributes[i], str) and level == k+1:
                dataframe.at[i,k+1] = meta_attributes[i]
    row_index_list = list(dataframe.count(axis=1) - 1)
    indices = [i for i, x in enumerate(row_index_list) if x == level]
    indices = list(dataframe.filter(items=indices, axis=0).index)
    level_params = []
    for index in indices:
        param = '/'.join(list(dataframe.loc[index].dropna()))
        if param.endswith(']'):
            level_params.append(param.rsplit('[',1)[0])
        else:
            level_params.append(param)

    return level_params


def get_params_dataframe(meta_params):
    """Get parameter dataframe"""

    level_count = 0
    row_count = len(meta_params)
    for i in range(row_count):
        meta_element = meta_params[i].split('/')
        level_count = \
            len(meta_element) if len(meta_element) > level_count else level_count

    dataframe = pd.DataFrame(index=range(row_count),
        columns=range(level_count-1), dtype=int)
    for i in range(row_count):
        meta_element = meta_params[i].split('/')
        for k in range(level_count):
            if k < len(meta_element):
                if meta_element[k].endswith(']'):
                    dataframe.at[i,k] = int(meta_element[k].split('[')[1][:-1])
                else:
                    dataframe.at[i,k] = 0
            else:
                dataframe.at[i,k] = None

    return dataframe


def get_dict_list_params(dict_list, old_level, new_level):
    """Get dictionary list parameters"""

    dict_parents = []
    dict_params = []
    for dict_element in dict_list:
        dict_parents.append(dict_element['path'] \
            .rsplit('/',old_level-new_level)[0])
        dict_params.append(dict_element['path'])

    return (dict_parents, dict_params)


def params2dict_list(dataframe, meta_params, meta_values, level):
    """Convert parameters to dictionary list"""

    dict_list = []
    row_index = list(dataframe.count(axis=1) - 1)
    col_index = list(dataframe.sum(axis=0))
    for i in range(len(col_index)):
        if col_index[i] > 0:
            row_index.append(i)
    row_index = pd.Series(row_index)
    dataframe_column = dataframe[level].dropna()
    row_level_index = row_index[row_index==level]
    dataframe_column = \
        dataframe_column.filter(items=list(row_level_index.index))
    dataframe_column_length = len(dataframe_column)
    if dataframe_column_length > 0:
        min_list = int(dataframe_column.min())
        max_list = int(dataframe_column.max())
        for list_index in range(min_list, max_list+1):
            dataframe_list = dataframe_column[dataframe_column == list_index]
            dataframe_index_list = \
                list(get_consecutive_list(dataframe_list.index.tolist()))
            for dataframe_indices in dataframe_index_list:
                for dataframe_index in dataframe_indices:
                    meta_element = meta_params[dataframe_index].split('/')
                    path = meta_params[dataframe_index]
                    parent = {}
                    parent['level'] = len(path.split('/'))-1
                    parent['path'] = path
                    parent['type'] = 'value'
                    parent['key'] = meta_element[level]
                    parent['list'] = list_index
                    parent['dict'] = meta_values[dataframe_index]
                    dict_list.append(parent)

    return dict_list


def upgrade_dictionary2level(old_dict_list, old_level, new_level):
    """Upgrade dictionary to level"""

    new_dict_list = []
    (dict_parents, _) = \
        get_dict_list_params(old_dict_list, old_level, new_level)
    parents = list(dict.fromkeys(dict_parents))

    for parent in parents:
        indices = [i for i, x in enumerate(dict_parents) if x == parent]
        for index in indices:
            parent_dict = {}
            level_diff = old_dict_list[index]['level'] - new_level
            path = old_dict_list[index]['path']
            old_dict = old_dict_list[index]['dict']
            for _ in range(level_diff):
                key = path.rsplit('/',1)[1]
                if ('Code' in key) and (key != 'postalCode'):
                    value = old_dict
                    old_dict = {}
                    if 'EOS' in key:
                        code_url = 'http://earthdata.nasa.gov/metadata/' \
                            'resources/Codelist.xml#' + key
                    else:
                        code_url = 'https://cdn.earthdata.nasa.gov/iso/' \
                            'resources/Codelist/gmxCodelists.xml#' + key
                    old_dict['@codeList'] = code_url
                    old_dict['@codeListValue'] = value
                    old_dict['value'] = value
                new_dict = {}
                new_dict[key.split('[')[0]] = old_dict
                path = path.rsplit('/',1)[0]
                old_dict = new_dict
            parent_dict['level'] = new_level
            if path.endswith(']'):
                list_index = int(path.rsplit('[',1)[1][:-1])
                path = path.rsplit('[',1)[0]
                parent_dict['path'] = path
                parent_dict['key'] = path.rsplit('/',1)[1]
                parent_dict['list'] = list_index
                parent_dict['type'] = 'list'
            else:
                if isinstance(new_dict, dict):
                    parent_dict['path'] = path
                elif isinstance(new_dict, (int, float, str)):
                    parent_dict['path'] = path.rsplit('/',1)[0]
                if '/' in path:
                    parent_dict['key'] = path.rsplit('/',1)[1]
                else:
                    parent_dict['key'] = path
                parent_dict['list'] = 0
                parent_dict['type'] = 'single'
                parent_dict['dict'] = new_dict
            if level_diff > 0:
                parent_dict['dict'] = new_dict
            else:
                parent_dict['dict'] = old_dict_list[index]['dict']
            new_dict_list.append(parent_dict)

    return new_dict_list


def merge_dictionary_params(old_dict_list, level_params, old_level, new_level):
    """Merge ditionary parameters"""

    merged_dict_list = []
    (_, dict_params) = \
        get_dict_list_params(old_dict_list, old_level, new_level)
    params = list(dict.fromkeys(dict_params))
    params = [i for i, group in groupby(level_params)]

    for param in params:
        indices = [i for i, x in enumerate(dict_params) if x == param]
        list_indices = []
        for index in indices:
            list_indices.append(old_dict_list[index]['list'])
        if len(list_indices) > 0:
            max_index = max(list_indices)
        else:
            max_index = 0

        ### list case
        if max_index > 0:
            dict_indices = [[] for i in range(max_index)]
            for index in indices:
                dict_indices[old_dict_list[index]['list']-1].append(index)
            num_indices = len(dict_indices[0])
            combine_dict = []
            for dict_index in dict_indices:
                old_dict = old_dict_list[dict_index[0]]
                if num_indices > 1:
                    list_dict = {}
                    for list_index in dict_index:
                        old_dict = old_dict_list[list_index]
                        key, value = next(iter(old_dict['dict'].items()))
                        list_dict[key] = value
                    combine_dict.append(list_dict)
                else:
                    combine_dict.append(old_dict['dict'])

        ### single value case
        else:
            combine_dict = {}
            dict_type = None
            for index in indices:
                old_dict = old_dict_list[index]
                if isinstance(old_dict['dict'], dict):
                    key, value = next(iter(old_dict['dict'].items()))
                    combine_dict[key] = value
                    dict_type = 'dict'
                elif isinstance(old_dict['dict'], (float, int, str)):
                    if dict_type == 'dict':
                        combine_dict['value'] = old_dict['dict']
                    else:
                        combine_dict = old_dict['dict']

        merged_dict = {}
        merged_dict['level'] = old_dict['level']
        if old_dict['key'] == 'type' and \
            'acquisitionInformation' in old_dict['path']:
            merged_dict['path'] = param
            merged_dict['key'] = old_dict['key']
        elif old_dict['key'] in ['id','href','nilReason','type','uom'] and \
            'additionalAttribute' not in old_dict['path']:
            merged_dict['path'] = rreplace(param, '/','/@')
            merged_dict['key'] = '@' + old_dict['key']
        else:
            merged_dict['path'] = param
            merged_dict['key'] = old_dict['key']
        merged_dict['list'] = old_dict['list']
        merged_dict['type'] = 'single'
        merged_dict['dict'] = combine_dict
        merged_dict_list.append(merged_dict)

    return merged_dict_list


def create_postgis_engine(config_file):
    """Create PostGIS engine"""

    config = configparser.ConfigParser()
    config.read(config_file)
    connection_string = \
        "postgresql://" + config.get('postgres', 'user') + \
        ":'" + config.get('postgres', 'pass') + "'" + \
        "@" + config.get('postgres', 'host') + ":5432" + \
        "/" + config.get('postgres', 'db')
    engine = create_engine(connection_string)

    return engine


def access_remote_zip(zip_file):
    """Access remote zip file"""

    ### Set up cookie jar
    login_url = 'https://urs.earthdata.nasa.gov/oauth/authorize?client_id=' \
        'BO_n7nTIlMljdvU6kRRB3g&response_type=code&redirect_uri=' \
        'https://auth.asf.alaska.edu/login&app_type=401'
    jar = requests.cookies.RequestsCookieJar()
    username = os.getenv('EDL_USER')
    password = os.getenv('EDL_PASS')
    login = requests.get(login_url, auth=(username, password))
    for row in login.history:
        jar.update(row.cookies.get_dict())

    zip_filehandle = remotezip.RemoteZip(zip_file, cookies=jar)

    return zip_filehandle
