#!/usr/bin/python3

"""Import packages"""
import argparse
from argparse import RawTextHelpFormatter
import sys
from metadata.util import access_remote_zip
from metadata.sentinel_metadata import get_sentinel_bursts


def sentinel_bursts2geojson(zip_url, burst_map_file, out_file):
    """Extract Sentinel burst metadata and save in a GeoJSON file"""

    ### Work out access to remote zip file
    zip_handle = access_remote_zip(zip_url)
    name_list = zip_handle.namelist()
    safe_dir = name_list[0]

    ### Getting Sentinel-1 burst information
    print('\nSentinel-1 bursts to GeoJSON ...')
    gdf = get_sentinel_bursts(safe_dir, burst_map_file, zip_handle, zip_url)
    zip_handle.close()
    gdf.to_file(out_file, driver = 'GeoJSON')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='sentinel_bursts2geojson',
        description='Extract burst information from a remote granule zip ' \
            'file and save it into GeoJSON file',
        formatter_class=RawTextHelpFormatter)
    parser.add_argument('zip_url', metavar='<zip file URL>',
        help='URL to the remote granule zip file')
    parser.add_argument('burst_map', metavar='<burst map file>',
        help='name of the burst map GeoJSON file')
    parser.add_argument('out_file', metavar='<output file>',
        help='name of the burst GeoJSON file')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()

    sentinel_bursts2geojson(args.zip_url, args.burst_map, args.out_file)
