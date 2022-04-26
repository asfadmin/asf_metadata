# Make sure GDAL is setup
try:
    from osgeo import ogr
except ImportError:
    raise ("""
    ERROR: Could not find the GDAL/OGR Python library bindings. 
    On Debian based systems you can install it with these commands:
        sudo apt install python-gdal libgdal-dev
        pip3 install gdal==$(gdal-config --version) --global-option=build_ext --global-option='-I/usr/include/gdal/'
    """)

# Import what the user is supposed to use directly
from asf_metadata import \
    generate_iso_metadata, \
    generate_iso_template, \
    iso_metadata_xml2json
