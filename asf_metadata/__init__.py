# to change color of output to user:
class colors:
    WARNING = '\033[93m'
    ERROR = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'

# Make sure GDAL is setup
try:
    from osgeo import ogr
except ImportError as e:
    raise ImportError(f"""
    {colors.ERROR}
    ERROR: Could not find the GDAL/OGR Python library bindings. 
    On Debian based systems you can install it with these commands:
        sudo apt install -y python3-gdal libgdal-dev
        python3 -m pip install gdal==$(gdal-config --version) --global-option=build_ext --global-option='-I/usr/include/gdal/'
    {colors.END}
    """) from e

# Import what the user is supposed to use directly
from asf_metadata.generate_iso_metadata import generate_iso_metadata
from asf_metadata.generate_iso_template import generate_iso_template
from asf_metadata.iso_metadata_xml2json import iso_metadata_xml2json
