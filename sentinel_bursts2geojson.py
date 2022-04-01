#!/usr/bin/python3

import argparse
from argparse import RawTextHelpFormatter
import sys
from metadata.util import accessRemoteZip
from metadata.sentinelMetadata import getSentinelBursts


def sentinel_bursts2geojson(zipUrl, burstMapFile, outFile):

  ### Work out access to remote zip file
  zf = accessRemoteZip(zipUrl)
  nameList = zf.namelist()
  safeDir = nameList[0]

  ### Getting Sentinel-1 burst information
  print('\nSentinel-1 bursts to GeoJSON ...')
  gdf = getSentinelBursts(safeDir, burstMapFile, zf, zipUrl)
  zf.close()
  gdf.to_file(outFile, driver = 'GeoJSON')


if __name__ == '__main__':

  parser = argparse.ArgumentParser(prog='sentinel_bursts2geojson',
    description='Extract burst information from a remote granule zip file ' \
      'and save it into GeoJSON file',
    formatter_class=RawTextHelpFormatter)
  parser.add_argument('zipUrl', metavar='<zip file URL>',
    help='URL to the remote granule zip file')
  parser.add_argument('burstMap', metavar='<burst map file>',
    help='name of the burst map GeoJSON file')
  parser.add_argument('outFile', metavar='<output file>',
    help='name of the burst GeoJSON file')
  if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
  args = parser.parse_args()

  sentinel_bursts2geojson(args.zipUrl, args.burstMap, args.outFile)
