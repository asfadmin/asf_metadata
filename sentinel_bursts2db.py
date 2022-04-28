#!/usr/bin/python3

"""Import packages"""
import argparse
from argparse import RawTextHelpFormatter
import sys
from metadata.util import create_postgis_engine, access_remote_zip
from metadata.sentinel_metadata import get_sentinel_bursts


def sentinel_bursts2db(zip_url, burst_map_file, config_file):
    """Extract Sentinel burst metadata and save in a database"""

    ### Connect to postreSQL database
    engine = create_postgis_engine(config_file)

    ### Work out access to remote zip file
    zip_handle = access_remote_zip(zip_url)
    name_list = zip_handle.namelist()
    safe_dir = name_list[0]

    ### Get Sentinel-1 file information
    print('\nSentinel-1 bursts to database ...')
    gdf = get_sentinel_bursts(safe_dir, burst_map_file, zip_handle, zip_url)
    zip_handle.close()
    gdf.to_postgis('sentinel_bursts', con=engine, if_exists='append',
      index=False)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='sentinel_bursts2db',
        description='Extract burst information from a remote granule zip ' \
            'file and save it into database',
        formatter_class=RawTextHelpFormatter)
    parser.add_argument('zip_url', metavar='<zip file URL>',
        help='URL to the remote granule zip file')
    parser.add_argument('burst_map', metavar='<burst map file>',
        help='name of the burst map GeoJSON file')
    parser.add_argument('config', metavar='<DB config file>',
        help='name of the database configuration file')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()

    sentinel_bursts2db(args.zip_url, args.burst_map, args.config)
