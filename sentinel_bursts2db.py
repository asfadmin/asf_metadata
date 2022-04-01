#!/usr/bin/python3

import argparse
from argparse import RawTextHelpFormatter
import sys
from metadata.util import create_postgis_engine, accessRemoteZip
from metadata.sentinelMetadata import getSentinelBursts


def sentinel_bursts2db(zipUrl, burstMapFile, configFile):

  ### Connect to postreSQL database
  engine = create_postgis_engine(configFile)

  ### Work out access to remote zip file
  zf = accessRemoteZip(zipUrl)
  nameList = zf.namelist()
  safeDir = nameList[0]

  ### Get Sentinel-1 file information
  print('\nSentinel-1 bursts to database ...')
  gdf = getSentinelBursts(safeDir, burstMapFile, zf, zipUrl)
  zf.close()
  gdf.to_postgis('sentinel_bursts', con=engine, if_exists='append', 
    index=False)


if __name__ == '__main__':

  parser = argparse.ArgumentParser(prog='sentinel_bursts2db',
    description='Extract burst information from a remote granule zip file ' \
      'and save it into database',
    formatter_class=RawTextHelpFormatter)
  parser.add_argument('zipUrl', metavar='<zip file URL>',
    help='URL to the remote granule zip file')
  parser.add_argument('burstMap', metavar='<burst map file>',
    help='name of the burst map GeoJSON file')
  parser.add_argument('config', metavar='<DB config file>',
    help='name of the database configuration file')
  if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
  args = parser.parse_args()

  sentinel_bursts2db(args.zipUrl, args.burstMap, args.config)
