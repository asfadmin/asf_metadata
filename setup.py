"""asf_metadata setuptools configuration"""
from setuptools import find_packages, setup

requirements = [
    "lxml",
    "numpy",
    "scipy",
    "openpyxl",
    "pandas",
    "geopandas",
    "sqlalchemy",
    "requests",
    "remotezip",
    "geoalchemy2",
    "psycopg2",
    "xmltodict",
]

test_requirements = [
    # pytest/etc goes here, once the test suite is setup.
]

with open("README.md", "r") as readme_file:
    readme = readme_file.read()

setup(
    name="asf_metadata",
    # version=Declared in pyproject.toml, through "[tool.setuptools_scm]"
    author="Alaska Satellite Facility Coop Team",
    author_email="UAF-ASF-CD@alaska.edu",
    description="Python wrapper for ASF's metadata",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/asfadmin/asf_metadata.git",
    project_urls={
        'Documentation': 'https://github.com/asfadmin/asf_metadata'
    },
    packages=find_packages(),
    package_dir={'asf_metadata': 'asf_metadata'},
    # Tells pypi to include these file types in the package
    package_data={"asf_metadata": ["*.png", "*.xlsx"]},
    python_requires='>=3.6',
    install_requires=requirements,
    extras_require={ "test": test_requirements },
    license='BSD',
    license_files=('LICENSE',),
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Topic :: Software Development",
        "Topic :: Scientific/Engineering :: Atmospheric Science",
        "Topic :: Scientific/Engineering :: GIS",
        "Topic :: Scientific/Engineering :: Hydrology",
        "Topic :: Utilities"
    ],
)
