#!/usr/bin/python3

"""Import packages"""
import os
from datetime import datetime, timedelta
import json
import lxml.etree as et
from osgeo import gdal, osr, ogr
import pandas as pd
import geopandas as gpd
import numpy as np
import scipy.interpolate
from shapely import wkt
import requests
import xmltodict


def xml2string(in_file, zip_handle=None):
    """Convert XML to string"""

    if zip_handle is None:
        parser = et.XMLParser(remove_blank_text=True)
        doc = et.parse(in_file, parser)
        string = et.tostring(doc)
    else:
        string = zip_handle.read(in_file)

    return string


def get_latlon_extent(filename):
    """Get lat/lon extent"""

    src = gdal.Open(filename)
    ulx, xres, _, uly, _, yres  = src.GetGeoTransform()
    lrx = ulx + (src.RasterXSize * xres)
    lry = uly + (src.RasterYSize * yres)

    source = osr.SpatialReference()
    source.ImportFromWkt(src.GetProjection())
    target = osr.SpatialReference()
    target.ImportFromEPSG(4326)
    transform = osr.CoordinateTransformation(source, target)

    (lon1, lat1, _) = transform.TransformPoint(ulx, uly)
    (lon2, lat2, _) = transform.TransformPoint(lrx, uly)
    (lon3, lat3, _) = transform.TransformPoint(ulx, lry)
    (lon4, lat4, _) = transform.TransformPoint(lrx, lry)

    lat_min = min(lat1,lat2,lat3,lat4)
    lat_max = max(lat1,lat2,lat3,lat4)
    lon_min = min(lon1,lon2,lon3,lon4)
    lon_max = max(lon1,lon2,lon3,lon4)

    return lat_min, lat_max, lon_min, lon_max


def extract_metadata(safe_dir, zip_handle=None):
    """Extract metadata"""

    manifest = read_manifest_file(os.path.join(safe_dir, 'manifest.safe'),
        zip_handle)
    sentinel_files = get_sentinel_files(safe_dir, manifest)
    file_count = len(sentinel_files)
    for i in range(file_count):
        file_name = sentinel_files[i]['name']
        file_type = sentinel_files[i]['type']
        if file_type == 'manifest':
            sentinel_files[i]['contents'] = \
                read_manifest_file(file_name, zip_handle)
        elif file_type == 'annotation':
            sentinel_files[i]['contents'] = \
                read_annotation_file(file_name, zip_handle)
        elif file_type == 'calibration':
            sentinel_files[i]['contents'] = \
                read_calibration_file(file_name, zip_handle)
        elif file_type == 'noise':
            sentinel_files[i]['contents'] = \
                read_noise_file(file_name, zip_handle)
        elif file_type == 'rfi':
            sentinel_files[i]['contents'] = \
                read_rfi_file(file_name, zip_handle)

    return sentinel_files


def extract_metadata_short(safe_dir, zip_handle=None):
    """Extract metadata"""

    manifest = read_manifest_file(os.path.join(safe_dir, 'manifest.safe'),
        zip_handle)
    sentinel_files = get_sentinel_files(safe_dir, manifest)
    file_count = len(sentinel_files)
    for i in range(file_count):
        file_name = sentinel_files[i]['name']
        file_type = sentinel_files[i]['type']
        if file_type == 'manifest':
            sentinel_files[i]['contents'] = \
                read_manifest_file(file_name, zip_handle)
        elif file_type == 'annotation':
            sentinel_files[i]['contents'] = \
                read_annotation_file(file_name, zip_handle)

    return sentinel_files


def read_manifest_file(manifest_file, zip_handle=None):
    """Read manifest file"""

    ns_xfdu = {'xfdu': 'urn:ccsds:schema:xfdu:1'}
    ns_all = dict(
        list(ns_xfdu.items())
    )

    ### Read manifest file
    doc = et.fromstring(xml2string(manifest_file, zip_handle))

    ### Pack first level into strings
    meta = {}
    meta['informationPackageMap'] = et.tostring( \
        doc.xpath('/xfdu:XFDU/informationPackageMap', namespaces=ns_all)[0])
    meta['metadataSection'] = \
        et.tostring(doc.xpath('/xfdu:XFDU/metadataSection',
            namespaces=ns_all)[0])
    meta['dataObjectSection'] = \
        et.tostring(doc.xpath('/xfdu:XFDU/dataObjectSection',
            namespaces=ns_all)[0])
    #print(meta['dataObjectSection'])

    return meta


def update_manifest_file(params):
    """Update manifest file"""

    meta = params

    return meta


def write_manifest_file(meta, manifest_file):
    """Write manifest file"""

    xfdu = '{urn:ccsds:schema:xfdu:1}'
    ns_xfdu = {'xfdu': 'urn:ccsds:schema:xfdu:1'}
    namespaces = dict(
        list(ns_xfdu.items())
    )

    # Build tree structure
    root = et.Element(f'{xfdu}XFDU', nsmap=namespaces)
    root.attrib['version'] = \
        'esa/safe/sentinel-1.0/sentinel-1/sar/level-1/slc/standard/iwdp'
    info_pack_map = et.fromstring(meta['informationPackageMap'])
    root.append(info_pack_map)
    meta_section = et.fromstring(meta['metadataSection'])
    root.append(meta_section)
    data_section = et.fromstring(meta['dataObjectSection'])
    root.append(data_section)

    # Write manifest file
    with open(manifest_file, 'wb') as file:
        file.write(et.tostring(root, xml_declaration=True, encoding='utf-8',
        pretty_print=True))


def read_annotation_file(annotation_file, zip_handle=None):
    """Read annotation file"""

    ### Read annotation file
    doc = et.fromstring(xml2string(annotation_file, zip_handle))

    ### Pack first level into strings
    meta = {}
    meta['adsHeader'] = et.tostring(doc.xpath('/product/adsHeader')[0])
    meta['qualityInformation'] = \
        et.tostring(doc.xpath('/product/qualityInformation')[0])
    meta['generalAnnotation'] = \
        et.tostring(doc.xpath('/product/generalAnnotation')[0])
    meta['imageAnnotation'] = \
        et.tostring(doc.xpath('/product/imageAnnotation')[0])
    meta['dopplerCentroid'] = \
        et.tostring(doc.xpath('/product/dopplerCentroid')[0])
    meta['antennaPattern'] = \
        et.tostring(doc.xpath('/product/antennaPattern')[0])
    meta['swathTiming'] = et.tostring(doc.xpath('/product/swathTiming')[0])
    meta['geolocationGrid'] = \
        et.tostring(doc.xpath('/product/geolocationGrid')[0])
    meta['coordinateConversion'] = \
        et.tostring(doc.xpath('/product/coordinateConversion')[0])
    meta['swathMerging'] = et.tostring(doc.xpath('/product/swathMerging')[0])

    return meta


def update_annotation(meta, burst_aoi, data_file):
    """Update nnnotation file"""

    return meta


def write_annotation_file(meta, annotation_file):
    """Write annotation file"""

    ### Put tree structure together
    root = et.Element('product')
    ads_header = et.fromstring(meta['adsHeader'])
    root.append(ads_header)
    quality_information = et.fromstring(meta['qualityInformation'])
    root.append(quality_information)
    general_annotation = et.fromstring(meta['generalAnnotation'])
    root.append(general_annotation)
    image_annotation = et.fromstring(meta['imageAnnotation'])
    root.append(image_annotation)
    doppler_centroid = et.fromstring(meta['dopplerCentroid'])
    root.append(doppler_centroid)
    antenna_pattern = et.fromstring(meta['antennaPattern'])
    root.append(antenna_pattern)
    swath_timing = et.fromstring(meta['swathTiming'])
    root.append(swath_timing)
    geolocation_grid = et.fromstring(meta['geolocationGrid'])
    root.append(geolocation_grid)
    coordinate_conversion = et.fromstring(meta['coordinateConversion'])
    root.append(coordinate_conversion)
    swath_merging = et.fromstring(meta['swathMerging'])
    root.append(swath_merging)

    # Write manifest file
    with open(annotation_file, 'wb') as file:
        file.write(et.tostring(root, xml_declaration=True, encoding='utf-8',
            pretty_print=True))


def annotation2spline(meta):
    """Extact spline function from annotation file"""

    ### Extract geolocation grid from metadata
    doc_grid = et.fromstring(meta['geolocationGrid'])
    number_of_points = int(doc_grid.xpath('/geolocationGrid/' \
        'geolocationGridPointList/@count')[0])

    ### Determine line/pixel indices for spline calculation
    lines = []
    pixels = []
    for point in range(number_of_points):
        xml = ('/geolocationGrid/geolocationGridPointList/' \
            f'geolocationGridPoint[{point+1}]/line')
        line = int(doc_grid.xpath(xml)[0].text)
        xml = ('/geolocationGrid/geolocationGridPointList/' \
            f'geolocationGridPoint[{point+1}]/pixel')
        pixel = int(doc_grid.xpath(xml)[0].text)
        if line == 0:
            pixels.append(pixel)
        if pixel == 0:
            lines.append(line)
    line_vector = np.array(lines)
    line_count = len(lines)
    pixel_vector = np.array(pixels)
    pixel_count = len(pixels)

    ### Build latitude/longitude grid
    lat_grid = np.zeros((line_count, pixel_count), dtype=float)
    lon_grid = np.zeros((line_count, pixel_count), dtype=float)
    height_grid = np.zeros((line_count, pixel_count), dtype=float)
    for point in range(number_of_points):
        xml = ('/geolocationGrid/geolocationGridPointList/' \
            f'geolocationGridPoint[{point+1}]/line')
        line = int(doc_grid.xpath(xml)[0].text)
        xml = ('/geolocationGrid/geolocationGridPointList/' \
            f'geolocationGridPoint[{point+1}]/pixel')
        pixel = int(doc_grid.xpath(xml)[0].text)
        xml = ('/geolocationGrid/geolocationGridPointList/' \
            f'geolocationGridPoint[{point+1}]/latitude')
        lat = float(doc_grid.xpath(xml)[0].text)
        xml = ('/geolocationGrid/geolocationGridPointList/' \
            f'geolocationGridPoint[{point+1}]/longitude')
        lon = float(doc_grid.xpath(xml)[0].text)
        xml = ('/geolocationGrid/geolocationGridPointList/' \
            f'geolocationGridPoint[{point+1}]/height')
        height = float(doc_grid.xpath(xml)[0].text)
        index_line = lines.index(line)
        index_pixel = pixels.index(pixel)
        lat_grid[index_line, index_pixel] = lat
        lon_grid[index_line, index_pixel] = lon
        height_grid[index_line, index_pixel] = height

    ### Calculate spline for interpolation
    lat_spline = scipy.interpolate.RectBivariateSpline(line_vector,
        pixel_vector, lat_grid)
    lon_spline = scipy.interpolate.RectBivariateSpline(line_vector,
        pixel_vector, lon_grid)
    height_spline = scipy.interpolate.RectBivariateSpline(line_vector,
        pixel_vector, height_grid)

    return (lat_spline, lon_spline, height_spline)


def annotation2burst_location(meta, burst_map, granule_info):
    """Extract burst location info from annoation file"""

    ### Initialize geodataframe and define data types
    dataframe = pd.DataFrame()

    (lat_spline, lon_spline, height_spline) = annotation2spline(meta)

    ### Extract metadata
    speed_of_light = 299792458.0
    doc = et.fromstring(meta['adsHeader'])
    swath = doc.xpath('/adsHeader/swath')[0].text
    swath_number = int(swath[-1:])
    polarization = doc.xpath('/adsHeader/polarisation')[0].text
    stop_time = doc.xpath('/adsHeader/stopTime')[0].text

    doc_general = et.fromstring(meta['generalAnnotation'])
    range_sampling_rate = float(doc_general.xpath('/generalAnnotation/' \
        'productInformation/rangeSamplingRate')[0].text)
    radar_wavelength = speed_of_light / \
        float(doc_general.xpath('/generalAnnotation/productInformation/' \
        'radarFrequency')[0].text)
    azimuth_steering_rate = float(doc_general.xpath('/generalAnnotation/' \
        'productInformation/azimuthSteeringRate')[0].text)
    prf = float(doc_general.xpath('/generalAnnotation/' \
        'downlinkInformationList/downlinkInformation/prf')[0].text)
    azimuth_fm_rate_list_count = int(doc_general.xpath('/generalAnnotation/' \
        'azimuthFmRateList/@count')[0])
    fm_rates = []
    for rate in range(azimuth_fm_rate_list_count):
        xml = ('/generalAnnotation/azimuthFmRateList/azimuthFmRate'
            f'[{rate+1}]/azimuthTime')
        ref_time = datetime.fromisoformat(doc_general.xpath(xml)[0].text)
        xml = ('/generalAnnotation/azimuthFmRateList/azimuthFmRate'
            f'[{rate+1}]/azimuthFmRatePolynomial')
        poly_coeffs = \
            list(map(float, doc_general.xpath(xml)[0].text.split(' ')))
        fm_rates.append((ref_time, poly_coeffs))

    doc_image = et.fromstring(meta['imageAnnotation'])
    azimuth_time_interval = float(doc_image.xpath('/imageAnnotation/' \
        'imageInformation/azimuthTimeInterval')[0].text)
    range_pixel_spacing = float(doc_image.xpath('/imageAnnotation/' \
        'imageInformation/rangePixelSpacing')[0].text)
    range_window_type = doc_image.xpath('/imageAnnotation/' \
        'processingInformation/swathProcParamsList/swathProcParams/' \
        'rangeProcessing/windowType')[0].text
    range_window_coefficient = float(doc_image.xpath('/imageAnnotation/' \
        'processingInformation/swathProcParamsList/swathProcParams/' \
        'rangeProcessing/windowCoefficient')[0].text)
    azimuth_window_type = doc_image.xpath('/imageAnnotation/' \
        'processingInformation/swathProcParamsList/swathProcParams/' \
        'azimuthProcessing/windowType')[0].text
    azimuth_window_coefficient = float(doc_image.xpath('/imageAnnotation/' \
        'processingInformation/swathProcParamsList/swathProcParams/' \
        'azimuthProcessing/windowCoefficient')[0].text)

    doc_doppler = et.fromstring(meta['dopplerCentroid'])
    dc_estimate_list_count = int(doc_doppler.xpath('/dopplerCentroid/' \
        'dcEstimateList/@count')[0])
    dopplers = []
    for doppler in range(1,dc_estimate_list_count+1):
        xml = (f'/dopplerCentroid/dcEstimateList/dcEstimate[{doppler}]/'
            'azimuthTime')
        ref_time = datetime.fromisoformat(doc_doppler.xpath(xml)[0].text)
        xml = (f'/dopplerCentroid/dcEstimateList/dcEstimate[{doppler}]/'
            'dataDcPolynomial')
        poly_coeffs = \
            list(map(float, doc_doppler.xpath(xml)[0].text.split(' ')))
        dopplers.append((ref_time, poly_coeffs))

    doc_swath = et.fromstring(meta['swathTiming'])
    lines_per_burst = \
        int(doc_swath.xpath('/swathTiming/linesPerBurst')[0].text)
    samples_per_burst = \
        int(doc_swath.xpath('/swathTiming/samplesPerBurst')[0].text)
    number_of_bursts = int(doc_swath.xpath('/swathTiming/burstList/@count')[0])

    doc_grid = et.fromstring(meta['geolocationGrid'])
    number_of_points = int(doc_grid.xpath('/geolocationGrid/' \
        'geolocationGridPointList/@count')[0])

    ### Determine line/pixel indices for spline calculation
    lines = []
    pixels = []
    slant_ranges = []
    for point in range(number_of_points):
        xml = ('/geolocationGrid/geolocationGridPointList/'
            f'geolocationGridPoint[{point+1}]/line')
        #print(doc_grid.xpath(xml))
        line = int(doc_grid.xpath(xml)[0].text)
        xml = ('/geolocationGrid/geolocationGridPointList/'
            f'geolocationGridPoint[{point+1}]/pixel')
        pixel = int(doc_grid.xpath(xml)[0].text)
        xml = ('/geolocationGrid/geolocationGridPointList/'
            f'geolocationGridPoint[{point+1}]/slantRangeTime')
        slant_range_time = float(doc_grid.xpath(xml)[0].text)
        if line == 0:
            pixels.append(pixel)
        if pixel == 0:
            lines.append(line)
        slant_range = slant_range_time * speed_of_light / 2.0
        slant_ranges.append(slant_range)
    max_pixel = max(pixels)
    center_pixel = max_pixel / 2.0

    acq_time = []
    sensing_time = []
    anx_time = []
    for burst in range(1,number_of_bursts+1):
        xml = (f'/swathTiming/burstList/burst[{burst}]/azimuthAnxTime')
        anx_time_burst = float(doc_swath.xpath(xml)[0].text)
        xml = (f'/swathTiming/burstList/burst[{burst}]/azimuthTime')
        azimuth_start_time = doc_swath.xpath(xml)[0].text
        xml = (f'/swathTiming/burstList/burst[{burst}]/sensingTime')
        sensing_start = doc_swath.xpath(xml)[0].text
        acq_time.append(datetime.fromisoformat(azimuth_start_time))
        sensing_time.append(datetime.fromisoformat(sensing_start))
        anx_time.append(anx_time_burst)
    acq_time.append(datetime.fromisoformat(stop_time))
    sensing_time.append(datetime.fromisoformat(stop_time))
    sensing_period = (lines_per_burst-1)/prf
    sensing_mid_period = sensing_period / 2.0

    ### Determine valid data and calculate burst boundaries
    for burst in range(number_of_bursts):

        first_line = lines[burst]
        last_line = lines[burst+1]
        starting_range = slant_ranges[burst]
        burst_time_period = \
            (last_line - first_line - 1) * azimuth_time_interval
        acq_period = abs(acq_time[burst+1] - acq_time[burst])
        acq_time_period = acq_period.seconds + acq_period.microseconds/1000000
        last_line += int(abs(burst_time_period - acq_time_period) / \
            azimuth_time_interval + 0.5)
        lats = np.zeros(4, dtype=float)
        lons = np.zeros(4, dtype=float)
        heights = np.zeros(4, dtype=float)
        lats[0] = lat_spline(first_line, 0)[0][0]
        lons[0] = lon_spline(first_line, 0)[0][0]
        heights[0] = height_spline(first_line, 0)[0][0]
        lats[1] = lat_spline(first_line, max_pixel)[0][0]
        lons[1] = lon_spline(first_line, max_pixel)[0][0]
        heights[1] = height_spline(first_line, max_pixel)[0][0]
        lats[2] = lat_spline(last_line, max_pixel)[0][0]
        lons[2] = lon_spline(last_line, max_pixel)[0][0]
        heights[2] = height_spline(last_line, max_pixel)[0][0]
        lats[3] = lat_spline(last_line, 0)[0][0]
        lons[3] = lon_spline(last_line, 0)[0][0]
        heights[3] = height_spline(last_line, 0)[0][0]
        center_line = first_line + (last_line - first_line) / 2.0
        terrain_height = height_spline(center_line, center_pixel)[0][0]

        xml = (f'/swathTiming/burstList/burst[{burst+1}]/azimuthAnxTime')
        anx_time_burst = float(doc_swath.xpath(xml)[0].text)
        xml = (f'/swathTiming/burstList/burst[{burst+1}]/azimuthTime')
        azimuth_start_time = doc_swath.xpath(xml)[0].text
        start_time = datetime.fromisoformat(azimuth_start_time)
        stop_time = start_time + timedelta(seconds=burst_time_period)
        azimuth_stop_time = stop_time.isoformat(timespec='microseconds')
        xml = (f'/swathTiming/burstList/burst[{burst+1}]/sensingTime')
        sensing_start = doc_swath.xpath(xml)[0].text
        start_time = datetime.fromisoformat(sensing_start)
        mid_time = start_time + timedelta(seconds=sensing_mid_period)
        stop_time = start_time + timedelta(seconds=sensing_period)
        sensing_stop = stop_time.isoformat(timespec='microseconds')

        time_diff = \
            [np.abs((mid_time - val[0]).total_seconds()) for val in fm_rates]
        index = np.argmin(time_diff)
        azimuth_fm_rate_polynomial = fm_rates[index][1]

        time_diff = \
            [np.abs((mid_time - val[0]).total_seconds()) for val in dopplers]
        index = np.argmin(time_diff)
        doppler_polynomial = dopplers[index][1]

        xml = (f'/swathTiming/burstList/burst[{burst+1}]/firstValidSample')
        valid_samples = \
            list(map(int, doc_swath.xpath(xml)[0].text.split(' ')))
        first_valid_sample = np.median(np.array(valid_samples)).astype(int)
        valid_indices = \
           [i for i, x in enumerate(valid_samples) if x == first_valid_sample]
        first_valid_line = int(min(valid_indices))
        last_valid_line = int(max(valid_indices))
        xml = (f'/swathTiming/burstList/burst[{burst+1}]/lastValidSample')
        valid_samples = list(map(int, doc_swath.xpath(xml)[0].text.split(' ')))
        last_valid_sample = np.median(np.array(valid_samples)).astype(int)
        number_valid_lines = last_valid_line - first_valid_line - 1
        number_valid_samples = last_valid_sample - first_valid_sample - 1

        ipf = float(granule_info['ipf'])
        if ipf >= 3.4:
            xml = (f'/swathTiming/burstList/burst[{burst+1}]/burstId/'
                '@absolute')
            absolute_esa_burst_id = int(doc_swath.xpath(xml)[0])
            xml = (f'/swathTiming/burstList/burst[{burst+1}]/burstId')
            relative_esa_burst_id = int(doc_swath.xpath(xml)[0].text)
        else:
            absolute_esa_burst_id = -1
            relative_esa_burst_id = -1

        orbit_direction = granule_info['passDirection']
        track = granule_info['trackNumber']
        absolute_burst_id = ('{get_burst_id(burst_map, orbit_direction, '
            'track, swath, float(anx_time[burst]))}-'
            f'{granule_info["orbitNumber"]:06d}')
        relative_burst_id = get_burst_id(burst_map, orbit_direction, track,
            swath, float(anx_time[burst]))
        if burst > 0:
            previous_burst_id = get_burst_id(burst_map, orbit_direction, track,
                swath, float(anx_time[burst-1]))
        else:
            previous_burst_id = None
        if burst < number_of_bursts - 1:
            next_burst_id = get_burst_id(burst_map, orbit_direction, track,
                swath, float(anx_time[burst+1]))
        else:
            next_burst_id = None

        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint_2D(float(lons[0]), float(lats[0]))
        ring.AddPoint_2D(float(lons[1]), float(lats[1]))
        ring.AddPoint_2D(float(lons[2]), float(lats[2]))
        ring.AddPoint_2D(float(lons[3]), float(lats[3]))
        ring.AddPoint_2D(float(lons[0]), float(lats[0]))
        polygon = ogr.Geometry(ogr.wkbPolygon)
        polygon.AddGeometry(ring)

        values = pd.DataFrame()
        values.at[0,'granuleUR'] = ("{granuleInfo['granule']}-BURSTS-{swath}"
            "{polarization}-{burst+1}:02d")
        values.at[0,'insertDate'] = None
        values.at[0,'productionDateTime'] = datetime.utcnow().isoformat() + 'Z'
        values.at[0,'identifier'] = ("{granuleInfo['granule']}-BURSTS")
        values.at[0,'inputGranule1'] = granule_info['granule']
        values.at[0,'beginningDateTime'] = azimuth_start_time
        values.at[0,'endingDateTime'] = azimuth_stop_time

        values.at[0,'platformShortName'] = granule_info['platform']
        values.at[0,'urlOfFrame'] = granule_info['urlOfFrame']
        values.at[0,'swath'] = swath
        values.at[0,'startLineInSwath'] = int(first_line)
        values.at[0,'endLineInSwath'] = int(last_line)
        values.at[0,'ipf'] = granule_info['ipf']

        values.at[0,'numberOfLines'] = int(lines_per_burst)
        values.at[0,'numberOfSamples'] = int(samples_per_burst)
        values.at[0,'startingRange'] = float(starting_range)
        values.at[0,'sensingStart'] = sensing_start
        values.at[0,'sensingStop'] = sensing_stop
        values.at[0,'burstStartUTC'] = azimuth_start_time
        values.at[0,'burstStopUTC'] = azimuth_stop_time
        values.at[0,'trackNumber'] = int(granule_info['trackNumber'])
        values.at[0,'frameNumber'] = int(-1)
        values.at[0,'orbitNumber'] = int(granule_info['orbitNumber'])
        values.at[0,'swathNumber'] = int(swath_number)
        values.at[0,'burstNumber'] = int(burst) + 1
        values.at[0,'passDirection'] = granule_info['passDirection']
        values.at[0,'azimuthSteeringRate'] = float(azimuth_steering_rate)
        values.at[0,'rangePixelSize'] = float(range_pixel_spacing)
        values.at[0,'rangeSamplingRate'] = float(range_sampling_rate)
        values.at[0,'azimuthTimeInterval'] = float(azimuth_time_interval)
        values.at[0,'radarWavelength'] = float(radar_wavelength)
        values.at[0,'polarization'] = polarization
        values.at[0,'terrainHeight'] = float(terrain_height)
        values.at[0,'prf'] = float(prf)
        values.at[0,'firstValidLine'] = int(first_valid_line)
        values.at[0,'numValidLines'] = int(number_valid_lines)
        values.at[0,'firstValidSample'] = int(first_valid_sample)
        values.at[0,'numValidSamples'] = int(number_valid_samples)
        values.at[0,'rangeWindowType'] = range_window_type
        values.at[0,'rangeWindowCoefficient'] = \
            float(range_window_coefficient)
        values.at[0,'azimuthWindowType'] = azimuth_window_type
        values.at[0,'azimuthWindowCoefficient'] = \
            float(azimuth_window_coefficient)
        values.at[0,'azimuthFMRate'] = \
            ','.join(map(str, azimuth_fm_rate_polynomial))
        values.at[0,'doppler'] = ','.join(map(str, doppler_polynomial))

        values.at[0,'absoluteBurstID'] = absolute_burst_id
        values.at[0,'relativeBurstID'] = relative_burst_id
        values.at[0,'previousBurstID'] = previous_burst_id
        values.at[0,'nextBurstID'] = next_burst_id
        values.at[0,'absoluteESAburstID'] = absolute_esa_burst_id
        values.at[0,'relativeESAburstID'] = relative_esa_burst_id
        values.at[0,'timeSinceAnxNode'] = anx_time_burst
        values.at[0,'ascendingNodeTime'] = granule_info['ascendingNodeTime']
        values.at[0,'hasLand'] = True
        values.at[0,'crossesDateLine'] = False
        values.at[0,'landFraction'] = 1.0

        values.at[0,'pointLongitude_1'] = float(lons[0])
        values.at[0,'pointLatitude_1'] = float(lats[0])
        values.at[0,'pointHeight_1'] = float(heights[0])
        values.at[0,'pointLongitude_2'] = float(lons[1])
        values.at[0,'pointLatitude_2'] = float(lats[1])
        values.at[0,'pointHeight_2'] = float(heights[1])
        values.at[0,'pointLongitude_3'] = float(lons[2])
        values.at[0,'pointLatitude_3'] = float(lats[2])
        values.at[0,'pointHeight_3'] = float(heights[2])
        values.at[0,'pointLongitude_4'] = float(lons[3])
        values.at[0,'pointLatitude_4'] = float(lats[3])
        values.at[0,'pointHeight_4'] = float(heights[3])
        values.at[0,'pointLongitude_5'] = float(lons[0])
        values.at[0,'pointLatitude_5'] = float(lats[0])
        values.at[0,'pointHeight_5'] = float(heights[0])
        values.at[0,'geometry'] = wkt.loads(polygon.ExportToWkt())

        dataframe = pd.concat([dataframe, values], ignore_index=True)

    return dataframe


def annotation2burst_location_short(meta, granule_info):
    """Extract burst location info from annoation file"""

    ### Initialize geodataframe and define data types
    dataframe = pd.DataFrame()

    (lat_spline, lon_spline, height_spline) = annotation2spline(meta)

    ### Extract metadata
    speed_of_light = 299792458.0
    doc = et.fromstring(meta['adsHeader'])
    swath = doc.xpath('/adsHeader/swath')[0].text
    stop_time = doc.xpath('/adsHeader/stopTime')[0].text

    doc_image = et.fromstring(meta['imageAnnotation'])
    azimuth_time_interval = float(doc_image.xpath('/imageAnnotation/' \
        'imageInformation/azimuthTimeInterval')[0].text)

    doc_swath = et.fromstring(meta['swathTiming'])
    number_of_bursts = int(doc_swath.xpath('/swathTiming/burstList/@count')[0])

    doc_grid = et.fromstring(meta['geolocationGrid'])
    number_of_points = int(doc_grid.xpath('/geolocationGrid/' \
        'geolocationGridPointList/@count')[0])

    ### Determine line/pixel indices for spline calculation
    lines = []
    pixels = []
    slant_ranges = []
    for point in range(number_of_points):
        xml = (f'/geolocationGrid/geolocationGridPointList/'
            f'geolocationGridPoint[{point+1}]/line')
        line = int(doc_grid.xpath(xml)[0].text)
        xml = (f'/geolocationGrid/geolocationGridPointList/'
            f'geolocationGridPoint[{point+1}]/pixel')
        pixel = int(doc_grid.xpath(xml)[0].text)
        xml = (f'/geolocationGrid/geolocationGridPointList/'
            f'geolocationGridPoint[{point+1}]/slantRangeTime')
        slant_range_time = float(doc_grid.xpath(xml)[0].text)
        if line == 0:
            pixels.append(pixel)
        if pixel == 0:
            lines.append(line)
        slant_range = slant_range_time * speed_of_light / 2.0
        slant_ranges.append(slant_range)
    max_pixel = max(pixels)

    acq_time = []
    anx_time = []
    for burst in range(1,number_of_bursts+1):
        xml = (f'/swathTiming/burstList/burst[{burst}]/azimuthAnxTime')
        anx_time_burst = float(doc_swath.xpath(xml)[0].text)
        xml = (f'/swathTiming/burstList/burst[{burst}]/azimuthTime')
        azimuth_start_time = doc_swath.xpath(xml)[0].text
        xml = (f'/swathTiming/burstList/burst[{burst}]/sensingTime')
        sensing_time = doc_swath.xpath(xml)[0].text
        acq_time.append(datetime.fromisoformat(azimuth_start_time))
        anx_time.append(anx_time_burst)
    acq_time.append(datetime.fromisoformat(stop_time))

    ### Determine valid data and calculate burst boundaries
    for burst in range(number_of_bursts):

        first_line = lines[burst]
        last_line = lines[burst+1]
        burst_time_period = \
            (last_line - first_line - 1) * azimuth_time_interval
        acq_period = abs(acq_time[burst+1] - acq_time[burst])
        acq_time_period = acq_period.seconds + acq_period.microseconds/1000000
        last_line += int(abs(burst_time_period - acq_time_period) / \
            azimuth_time_interval + 0.5)
        lats = np.zeros(4, dtype=float)
        lons = np.zeros(4, dtype=float)
        heights = np.zeros(4, dtype=float)
        lats[0] = lat_spline(first_line, 0)[0][0]
        lons[0] = lon_spline(first_line, 0)[0][0]
        heights[0] = height_spline(first_line, 0)[0][0]
        lats[1] = lat_spline(first_line, max_pixel)[0][0]
        lons[1] = lon_spline(first_line, max_pixel)[0][0]
        heights[1] = height_spline(first_line, max_pixel)[0][0]
        lats[2] = lat_spline(last_line, max_pixel)[0][0]
        lons[2] = lon_spline(last_line, max_pixel)[0][0]
        heights[2] = height_spline(last_line, max_pixel)[0][0]
        lats[3] = lat_spline(last_line, 0)[0][0]
        lons[3] = lon_spline(last_line, 0)[0][0]
        heights[3] = height_spline(last_line, 0)[0][0]

        xml = (f'/swathTiming/burstList/burst[{burst+1}]/azimuthAnxTime')
        anx_time_burst = float(doc_swath.xpath(xml)[0].text)
        xml = (f'/swathTiming/burstList/burst[{burst+1}]/azimuthTime')
        azimuth_start_time = doc_swath.xpath(xml)[0].text
        xml = (f'/swathTiming/burstList/burst[{burst+1}]/sensingTime')
        sensing_time = doc_swath.xpath(xml)[0].text
        start_time = datetime.fromisoformat(azimuth_start_time)
        stop_time = start_time + timedelta(seconds=burst_time_period)
        azimuth_stop_time = stop_time.isoformat(timespec='microseconds')

        ipf = float(granule_info['ipf'])
        if ipf >= 3.4:
            xml = (f'/swathTiming/burstList/burst[{burst+1}]/burstId/'
                '@absolute')
            absolute_esa_burst_id = int(doc_swath.xpath(xml)[0])
            xml = (f'/swathTiming/burstList/burst[{burst+1}]/burstId')
            relative_esa_burst_id = int(doc_swath.xpath(xml)[0].text)
        else:
            absolute_esa_burst_id = -1
            relative_esa_burst_id = -1

        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint_2D(float(lons[0]), float(lats[0]))
        ring.AddPoint_2D(float(lons[1]), float(lats[1]))
        ring.AddPoint_2D(float(lons[2]), float(lats[2]))
        ring.AddPoint_2D(float(lons[3]), float(lats[3]))
        ring.AddPoint_2D(float(lons[0]), float(lats[0]))
        polygon = ogr.Geometry(ogr.wkbPolygon)
        polygon.AddGeometry(ring)

        values = pd.DataFrame()
        values.at[0,'granule'] = granule_info['granule']
        values.at[0,'swath'] = swath
        values.at[0,'trackNumber'] = int(granule_info['trackNumber'])
        values.at[0,'orbitNumber'] = int(granule_info['orbitNumber'])
        values.at[0,'burstNumber'] = int(burst) + 1
        values.at[0,'orbitDirection'] = granule_info['passDirection']
        values.at[0,'absoluteESAburstID'] = absolute_esa_burst_id
        values.at[0,'relativeESAburstID'] = relative_esa_burst_id
        values.at[0,'burstStartUTC'] = azimuth_start_time
        values.at[0,'burstStopUTC'] = azimuth_stop_time
        values.at[0,'burstSensingTime'] = sensing_time
        values.at[0,'startLineInSwath'] = int(first_line)
        values.at[0,'endLineInSwath'] = int(last_line)
        values.at[0,'ipf'] = granule_info['ipf']
        values.at[0,'timeSinceAnxNode'] = anx_time_burst
        values.at[0,'ascendingNodeTime'] = granule_info['ascendingNodeTime']
        values.at[0,'geometry'] = wkt.loads(polygon.ExportToWkt())

        dataframe = pd.concat([dataframe, values], ignore_index=True)

    return dataframe


def read_calibration_file(calibration_file, zf=None):
    """Read calibration file"""

    # Read annotation file
    doc = et.fromstring(xml2string(calibration_file, zf))

    ### Pack first level into strings
    meta = {}
    meta['adsHeader'] = et.tostring(doc.xpath('/calibration/adsHeader')[0])
    meta['calibrationInformation'] = \
        et.tostring(doc.xpath('/calibration/calibrationInformation')[0])
    meta['calibrationVectorList'] = \
        et.tostring(doc.xpath('/calibration/calibrationVectorList')[0])

    return meta


def update_calibration(meta, burst_aoi):
    """Update calibration information"""

    return meta


def write_calibration_file(meta, calibration_file):
    """Write calibration file"""

    # Put tree structure together
    root = et.Element('calibration')
    ads_header = et.fromstring(meta['adsHeader'])
    root.append(ads_header)
    calibration_information = et.fromstring(meta['calibrationInformation'])
    root.append(calibration_information)
    calibration_vector_list = et.fromstring(meta['calibrationVectorList'])
    root.append(calibration_vector_list)

    # Write manifest file
    with open(calibration_file, 'wb') as file:
        file.write(et.tostring(root, xml_declaration=True, encoding='utf-8',
        pretty_print=True))

def read_noise_file(noise_file, zip_file=None):
    """Read noise file"""

    ### Read noise file
    doc = et.fromstring(xml2string(noise_file, zip_file))

    ### Pack first level into strings
    meta = {}
    meta['adsHeader'] = et.tostring(doc.xpath('/noise/adsHeader')[0])
    noise_vector_list = doc.xpath('/noise/noiseVectorList')
    noise_range_vector_list = doc.xpath('/noise/noiseRangeVectorList')
    noise_azimuth_vector_list = doc.xpath('/noise/noiseAzimuthVectorList')
    if len(noise_vector_list) > 0:
        meta['noiseVectorList'] = et.tostring(noise_vector_list[0])
    else:
        meta['noiseRangeVectorList'] = et.tostring(noise_range_vector_list[0])
        meta['noiseAzimuthVectorList'] = \
            et.tostring(noise_azimuth_vector_list[0])

    return meta


def update_noise(meta, burst_aoi):
    """Update noise information"""

    return meta


def write_noise_file(meta, noise_file):
    """Write noise file"""

    ### Put tree structure together
    root = et.Element('noise')
    ads_header = et.fromstring(meta['adsHeader'])
    root.append(ads_header)
    if 'noiseVectorList' in meta:
        noise_vector_list = et.fromstring(meta['noiseVectorList'])
        root.append(noise_vector_list)
    else:
        noise_range_vector_list = et.fromstring(meta['noiseRangeVectorList'])
        root.append(noise_range_vector_list)
        noise_azimuth_vector_list = \
            et.fromstring(meta['noiseAzimuthVectorList'])
        root.append(noise_azimuth_vector_list)

    ### Write noise file
    with open(noise_file, 'wb') as file:
        file.write(et.tostring(root, xml_declaration=True, encoding='utf-8',
        pretty_print=True))


def read_rfi_file(rfi_file, zip_file=None):
    """Read RFI file"""

    ### Read RFI file
    doc = et.fromstring(xml2string(rfi_file, zip_file))

    ### Pack first level into strings
    meta = {}
    meta['adsHeader'] = et.tostring(doc.xpath('/rfi/adsHeader')[0])
    meta['rfiMitigationApplied'] = \
        et.tostring(doc.xpath('/rfi/rfiMitigationApplied')[0])
    meta['rfiDetectionFromNoiseReportList'] = \
        et.tostring(doc.xpath('/rfi/rfiDetectionFromNoiseReportList')[0])
    meta['rfiBurstReportList'] = \
        et.tostring(doc.xpath('/rfi/rfiBurstReportList')[0])

    return meta


def update_rfi(meta, burst_aoi):
    """Update RFI information"""

    return meta


def write_rfi_file(meta, rfi_file):
    """Write RFI file"""

    ### Put tree structure together
    root = et.Element('rfi')
    ads_header = et.fromstring(meta['adsHeader'])
    root.append(ads_header)
    rfi_mitigation_applied = et.fromstring(meta['rfiMitigationApplied'])
    root.append(rfi_mitigation_applied)
    rfi_detection_from_noise_report_list = \
        et.fromstring(meta['rfiDetectionFromNoiseReportList'])
    root.append(rfi_detection_from_noise_report_list)
    rfi_burst_report_list = et.fromstring(meta['rfiBurstReportList'])
    root.append(rfi_burst_report_list)

    ### Write noise file
    with open(rfi_file, 'wb') as file:
        file.write(et.tostring(root, xml_declaration=True, encoding='utf-8',
        pretty_print=True))


def get_burst_id(burst_map, orbit_direction, track, swath, anx_time):
    """Get burst index"""

    burst_direction = burst_map[burst_map['orbitDir'] == orbit_direction]
    burst_track = burst_direction[burst_direction['track'] == track]
    burst_swath = pd.DataFrame(burst_track[burst_track['swath'] == swath])
    burst_swath = burst_swath.reset_index()
    delta_anx_time = 12*24*3600/175
    burst_id = 'TTT-AAAA-IIW-PP-VU'
    for i in range(len(burst_swath)):
        diff = abs(anx_time - burst_swath.at[i,'anxTime'])
        if diff < delta_anx_time:
            burst_id = burst_swath.at[i,'burstID']
            delta_anx_time = diff

    return burst_id


def get_sentinel_bursts(safe_dir, burst_map_file, zip_handle=None, url_str=None):
    """Get Sentinel bursts"""

    ### Read burst map
    base_name = os.path.basename(burst_map_file)
    message = f'Reading burst map file ({base_name}) ...'
    print(message)
    burst_map = gpd.read_file(burst_map_file)

    ### Define namespaces for XML parsing
    ns_safe = {'safe': 'http://www.esa.int/safe/sentinel-1.0'}
    ns_s1 = {'s1': 'http://www.esa.int/safe/sentinel-1.0/sentinel-1'}
    names = dict(
        list(ns_safe.items()) +
        list(ns_s1.items())
    )

    ### Extract granule metadata
    granule = safe_dir[:-6]
    print(f'Extracting information for granule ({granule}) ...')
    sentinel_files = extract_metadata(safe_dir, zip_handle)
    for sentinel_file in sentinel_files:
        if sentinel_file['type'] == 'manifest':
            manifest = sentinel_file['contents']
        if sentinel_file['type'] == 'annotation':
            if sentinel_file['polarization'] in ('VV','HH'):
                if sentinel_file['swath'] == 'IW1':
                    annotation_iw1 = sentinel_file['contents']
                elif sentinel_file['swath'] == 'IW2':
                    annotation_iw2 = sentinel_file['contents']
                elif sentinel_file['swath'] == 'IW3':
                    annotation_iw3 = sentinel_file['contents']
    '''
    (manifest, annotation_iw1, annotation_iw2, annotation_iw3) = \
        extract_metadata(safe_dir, zip_handle)
    '''

    doc = et.fromstring(manifest['metadataSection'])
    meta_object_count = len(doc.xpath('/metadataSection/metadataObject/@ID'))
    for i in range(1,meta_object_count+1):
        index = doc.xpath(f'/metadataSection/metadataObject[{i}]/@ID')[0]
        if index == 'processing':
            xml = (f'/metadataSection/metadataObject[{i}]/'
                'metadataWrap/xmlData/safe:processing/safe:facility/'
                'safe:software/@version')
            ipf = doc.xpath(xml, namespaces=names)[0]
        elif index == 'platform':
            xml = (f'/metadataSection/metadataObject[{i}]/'
                'metadataWrap/xmlData/safe:platform/safe:familyName')
            param = doc.xpath(xml, namespaces=names)
            platform = param[0].text
            xml = (f'/metadataSection/metadataObject[{i}]/'
                'metadataWrap/xmlData/safe:platform/safe:number')
            param = doc.xpath(xml, namespaces=names)
            platform += param[0].text
        elif index == 'measurementOrbitReference':
            xml = (f'/metadataSection/metadataObject[{i}]/'
                'metadataWrap/xmlData/safe:orbitReference/'
                'safe:orbitNumber[@type="start"]')
            param = doc.xpath(xml, namespaces=names)
            abolute_orbit_number = int(param[0].text)
            xml = (f'/metadataSection/metadataObject[{i}]/'
                'metadataWrap/xmlData/safe:orbitReference/'
                'safe:relativeOrbitNumber[@type="start"]')
            param = doc.xpath(xml, namespaces=names)
            relative_orbit_number = int(param[0].text)
            xml = (f'/metadataSection/metadataObject[{i}]/'
                'metadataWrap/xmlData/safe:orbitReference/safe:extension/'
                's1:orbitProperties/s1:pass')
            param = doc.xpath(xml, namespaces=names)
            orbit_direction = param[0].text
            xml = (f'/metadataSection/metadataObject[{i}]/'
                'metadataWrap/xmlData/safe:orbitReference/safe:extension/'
                's1:orbitProperties/s1:ascendingNodeTime')
            param = doc.xpath(xml, namespaces=names)
            ascending_node_time = param[0].text

    if not url_str:
        ### Determine URL (if needed)
        url = 'https://api.daac.asf.alaska.edu/services/search/' \
            'param?granule_list='
        url += safe_dir[:-6]
        data = requests.get(url).content
        data = json.dumps(xmltodict.parse(data.decode()))
        string = json.loads(data)
        url_str = string['metalink']['files']['file'][1]['resources']['url'] \
            ['#text']

    ### Split data files into bursts
    granule_info = {}
    granule_info['platform'] = platform
    granule_info['orbitNumber'] = abolute_orbit_number
    granule_info['relativeOrbitNumber'] = relative_orbit_number
    granule_info['trackNumber'] = relative_orbit_number
    granule_info['passDirection'] = orbit_direction
    granule_info['ascendingNodeTime'] = ascending_node_time
    granule_info['granule'] = safe_dir[:-6]
    granule_info['urlOfFrame'] = url_str
    granule_info['ipf'] = ipf

    #dataframe = initializeGeoDataFrame()
    dataframe = pd.DataFrame()

    ### Extraction burst information for all IW swaths
    burst_info = \
        annotation2burst_location(annotation_iw1, burst_map, granule_info)
    dataframe = pd.concat([dataframe, burst_info], ignore_index=True)
    burst_info = \
        annotation2burst_location(annotation_iw2, burst_map, granule_info)
    dataframe = pd.concat([dataframe, burst_info], ignore_index=True)
    burst_info = \
        annotation2burst_location(annotation_iw3, burst_map, granule_info)
    dataframe = pd.concat([dataframe, burst_info], ignore_index=True)
    dataframe = fix_data_types(dataframe)
    gdf = gpd.GeoDataFrame(dataframe, crs='EPSG:4326')
    gdf.set_geometry('geometry')

    return gdf


def get_sentinel_burst_info(safe_dir, zip_handle=None, url_str=None):
    """Get Sentinel bursts"""

    ### Define namespaces for XML parsing
    ns_safe = {'safe': 'http://www.esa.int/safe/sentinel-1.0'}
    ns_s1 = {'s1': 'http://www.esa.int/safe/sentinel-1.0/sentinel-1'}
    names = dict(
        list(ns_safe.items()) +
        list(ns_s1.items())
    )

    ### Extract granule metadata
    granule = safe_dir[:-6]
    print(f'Extracting information for granule ({granule}) ...')
    sentinel_files = extract_metadata_short(safe_dir, zip_handle)
    for sentinel_file in sentinel_files:
        if sentinel_file['type'] == 'manifest':
            manifest = sentinel_file['contents']
        if sentinel_file['type'] == 'annotation':
            if sentinel_file['polarization'] in ('VV','HH'):
                if sentinel_file['swath'] == 'IW1':
                    annotation_iw1 = sentinel_file['contents']
                elif sentinel_file['swath'] == 'IW2':
                    annotation_iw2 = sentinel_file['contents']
                elif sentinel_file['swath'] == 'IW3':
                    annotation_iw3 = sentinel_file['contents']
    '''
    (manifest, annotation_iw1, annotation_iw2, annotation_iw3) = \
        extract_metadata(safe_dir, zip_handle)
    '''

    doc = et.fromstring(manifest['metadataSection'])
    meta_object_count = len(doc.xpath('/metadataSection/metadataObject/@ID'))
    for i in range(1,meta_object_count+1):
        index = doc.xpath(f'/metadataSection/metadataObject[{i}]/@ID')[0]
        if index == 'processing':
            xml = (f'/metadataSection/metadataObject[{i}]/'
                'metadataWrap/xmlData/safe:processing/safe:facility/'
                'safe:software/@version')
            ipf = doc.xpath(xml, namespaces=names)[0]
        elif index == 'measurementOrbitReference':
            xml = (f'/metadataSection/metadataObject[{i}]/'
                'metadataWrap/xmlData/safe:orbitReference/'
                'safe:orbitNumber[@type="start"]')
            param = doc.xpath(xml, namespaces=names)
            absolute_orbit_number = int(param[0].text)
            xml = (f'/metadataSection/metadataObject[{i}]/'
                'metadataWrap/xmlData/safe:orbitReference/'
                'safe:relativeOrbitNumber[@type="start"]')
            param = doc.xpath(xml, namespaces=names)
            relative_orbit_number = int(param[0].text)
            xml = (f'/metadataSection/metadataObject[{i}]/'
                'metadataWrap/xmlData/safe:orbitReference/safe:extension/'
                's1:orbitProperties/s1:pass')
            param = doc.xpath(xml, namespaces=names)
            orbit_direction = param[0].text
            xml = (f'/metadataSection/metadataObject[{i}]/'
                'metadataWrap/xmlData/safe:orbitReference/safe:extension/'
                's1:orbitProperties/s1:ascendingNodeTime')
            param = doc.xpath(xml, namespaces=names)
            ascending_node_time = param[0].text

    if not url_str:
        ### Determine URL (if needed)
        url = 'https://api.daac.asf.alaska.edu/services/search/' \
            'param?granule_list='
        url += safe_dir[:-6]
        data = requests.get(url).content
        data = json.dumps(xmltodict.parse(data.decode()))
        string = json.loads(data)
        url_str = string['metalink']['files']['file'][1]['resources']['url'] \
            ['#text']

    ### Split data files into bursts
    granule_info = {}
    granule_info['orbitNumber'] = absolute_orbit_number
    granule_info['relativeOrbitNumber'] = relative_orbit_number
    granule_info['trackNumber'] = relative_orbit_number
    granule_info['passDirection'] = orbit_direction
    granule_info['ascendingNodeTime'] = ascending_node_time
    granule_info['granule'] = safe_dir[:-6]
    granule_info['urlOfFrame'] = url_str
    granule_info['ipf'] = ipf

    #dataframe = initializeGeoDataFrame()
    dataframe = pd.DataFrame()

    ### Extraction burst information for all IW swaths
    burst_info = \
        annotation2burst_location_short(annotation_iw1, granule_info)
    dataframe = pd.concat([dataframe, burst_info], ignore_index=True)
    burst_info = \
        annotation2burst_location_short(annotation_iw2, granule_info)
    dataframe = pd.concat([dataframe, burst_info], ignore_index=True)
    burst_info = \
        annotation2burst_location_short(annotation_iw3, granule_info)
    dataframe = pd.concat([dataframe, burst_info], ignore_index=True)
    dataframe = fix_data_types_short(dataframe)
    gdf = gpd.GeoDataFrame(dataframe, crs='EPSG:4326')
    gdf.set_geometry('geometry')

    return gdf


def fix_data_types(dataframe):
    """Fix data types"""

    dataframe['insertDate'] = pd.to_datetime(dataframe['insertDate'])
    dataframe['productionDateTime'] = \
        pd.to_datetime(dataframe['productionDateTime'])
    dataframe['beginningDateTime'] = \
        pd.to_datetime(dataframe['beginningDateTime'])
    dataframe['endingDateTime'] = pd.to_datetime(dataframe['endingDateTime'])

    dataframe['startLineInSwath'] = \
        dataframe['startLineInSwath'].astype('int32')
    dataframe['endLineInSwath'] = dataframe['endLineInSwath'].astype('int32')

    dataframe['numberOfLines'] = dataframe['numberOfLines'].astype('int32')
    dataframe['numberOfSamples'] = dataframe['numberOfSamples'].astype('int32')
    dataframe['startingRange'] = dataframe['startingRange'].astype('float32')
    dataframe['sensingStart'] = pd.to_datetime(dataframe['sensingStart'])
    dataframe['sensingStop'] = pd.to_datetime(dataframe['sensingStop'])
    dataframe['burstStartUTC'] = pd.to_datetime(dataframe['burstStartUTC'])
    dataframe['burstStopUTC'] = pd.to_datetime(dataframe['burstStopUTC'])
    dataframe['burstSensingTime'] = \
        pd.to_datetime(dataframe['burstSensingTime'])
    dataframe['trackNumber'] = dataframe['trackNumber'].astype('int32')
    dataframe['frameNumber'] = dataframe['frameNumber'].astype('int32')
    dataframe['orbitNumber'] = dataframe['orbitNumber'].astype('int32')
    dataframe['swathNumber'] = dataframe['swathNumber'].astype('int32')
    dataframe['burstNumber'] = dataframe['burstNumber'].astype('int32')
    dataframe['azimuthSteeringRate'] = \
        dataframe['azimuthSteeringRate'].astype('float32')
    dataframe['rangePixelSize'] = dataframe['rangePixelSize'].astype('float32')
    dataframe['rangeSamplingRate'] = \
        dataframe['rangeSamplingRate'].astype('float32')
    dataframe['azimuthTimeInterval'] = \
        dataframe['azimuthTimeInterval'].astype('float32')
    dataframe['radarWavelength'] = \
        dataframe['radarWavelength'].astype('float32')
    dataframe['terrainHeight'] = dataframe['terrainHeight'].astype('float32')
    dataframe['prf'] = dataframe['prf'].astype('float32')
    dataframe['firstValidLine'] = \
        dataframe['firstValidLine'].astype('int32')
    dataframe['numValidLines'] = dataframe['numValidLines'].astype('int32')
    dataframe['firstValidSample'] = \
        dataframe['firstValidSample'].astype('int32')
    dataframe['numValidSamples'] = \
        dataframe['numValidSamples'].astype('int32')
    dataframe['rangeWindowCoefficient'] = \
        dataframe['rangeWindowCoefficient'].astype('float32')
    dataframe['azimuthWindowCoefficient'] = \
        dataframe['azimuthWindowCoefficient'].astype('float32')

    dataframe['absoluteESAburstID'] = \
        dataframe['absoluteESAburstID'].astype('int32')
    dataframe['relativeESAburstID'] = \
        dataframe['relativeESAburstID'].astype('int32')
    dataframe['timeSinceAnxNode'] = \
        dataframe['timeSinceAnxNode'].astype('float32')
    dataframe['ascendingNodeTime'] = \
        pd.to_datetime(dataframe['ascendingNodeTime'])
    dataframe['hasLand'] = dataframe['hasLand'].astype('bool')
    dataframe['crossesDateLine'] = dataframe['crossesDateLine'].astype('bool')
    dataframe['landFraction'] = dataframe['landFraction'].astype('float32')

    dataframe['pointLongitude_1'] = \
        dataframe['pointLongitude_1'].astype('float32')
    dataframe['pointLatitude_1'] = \
        dataframe['pointLatitude_1'].astype('float32')
    dataframe['pointLongitude_2'] = \
        dataframe['pointLongitude_2'].astype('float32')
    dataframe['pointLatitude_2'] = \
        dataframe['pointLatitude_2'].astype('float32')
    dataframe['pointLongitude_3'] = \
        dataframe['pointLongitude_3'].astype('float32')
    dataframe['pointLatitude_3'] = \
        dataframe['pointLatitude_3'].astype('float32')
    dataframe['pointLongitude_4'] = \
        dataframe['pointLongitude_4'].astype('float32')
    dataframe['pointLatitude_4'] = \
        dataframe['pointLatitude_4'].astype('float32')
    dataframe['pointLongitude_5'] = \
        dataframe['pointLongitude_5'].astype('float32')
    dataframe['pointLatitude_5'] = \
        dataframe['pointLatitude_5'].astype('float32')

    return dataframe


def fix_data_types_short(dataframe):
    """Fix data types"""

    dataframe['trackNumber'] = dataframe['trackNumber'].astype('int32')
    dataframe['orbitNumber'] = dataframe['orbitNumber'].astype('int32')
    dataframe['burstNumber'] = dataframe['burstNumber'].astype('int32')
    dataframe['absoluteESAburstID'] = \
        dataframe['absoluteESAburstID'].astype('int32')
    dataframe['relativeESAburstID'] = \
        dataframe['relativeESAburstID'].astype('int32')
    dataframe['burstStartUTC'] = pd.to_datetime(dataframe['burstStartUTC'])
    dataframe['burstStopUTC'] = pd.to_datetime(dataframe['burstStopUTC'])
    dataframe['burstSensingTime'] = \
        pd.to_datetime(dataframe['burstSensingTime'])
    dataframe['startLineInSwath'] = \
        dataframe['startLineInSwath'].astype('int32')
    dataframe['endLineInSwath'] = dataframe['endLineInSwath'].astype('int32')
    dataframe['timeSinceAnxNode'] = \
        dataframe['timeSinceAnxNode'].astype('float32')
    dataframe['ascendingNodeTime'] = \
        pd.to_datetime(dataframe['ascendingNodeTime'])

    return dataframe


def parse_line(line, output_type):
    """Parse a line"""

    return_value = None
    if line.count(' - ') == 2:
        (datestamp, log_level, remainder) = line.split(' - ')
        if output_type == 'dateString':
            return_value = datetime.strptime(datestamp,
                '%m/%d/%Y %H:%M:%S %p').isoformat() + 'Z'
        elif output_type == 'log_level':
            return_value = log_level.split()
        if ':' in remainder:
            (key, value) = remainder.split(':', 1)
            if output_type == 'key':
                return_value = key.strip()
            elif output_type == 'value':
                return_value = value.strip()
    elif line.count(' - ') == 0 and line.count(': ') == 1:
        (key, value) = line.split(': ')
        if output_type == 'key':
            return_value = key
        elif output_type == 'value':
            return_value = value.strip()

    return return_value


def get_sentinel_files(safe, meta):
    """Get Sentinel file information"""

    files = []
    file = {}
    file['name'] = os.path.join(safe, 'manifest.safe')
    file['type'] = 'manifest'
    files.append(file)
    doc = et.fromstring(meta['dataObjectSection'])
    data_objects = doc.xpath('/dataObjectSection/dataObject')
    data_object_count = len(data_objects)
    for i in range(data_object_count):
        file = {}
        xml = (f'/dataObjectSection/dataObject[{i+1}]/byteStream/' \
            'fileLocation/@href')
        href = doc.xpath(xml)[0]
        file_parts = os.path.basename(href).split('-')
        file['name'] = os.path.join(safe, href[2:])
        xml = (f'/dataObjectSection/dataObject[{i+1}]/@ID')
        file['index'] = doc.xpath(xml)[0]
        if file['index'] == 'mapoverlay':
            file['type'] = 'preview'
        elif file['index'] == 'productpreview':
            file['type'] = 'preview'
        elif file['index'] == 'quicklook':
            file['type'] = 'preview'
        elif file['index'].startswith('product'):
            file['type'] = 'annotation'
            file['swath'] = file_parts[1].upper()
            file['polarization'] = file_parts[3].upper()
            file['start_time'] = file_parts[4]
            file['stop_time'] = file_parts[5]
        elif file['index'].startswith('calibration'):
            file['type'] = 'calibration'
            file['swath'] = file_parts[2].upper()
            file['polarization'] = file_parts[4].upper()
            file['start_time'] = file_parts[5]
            file['stop_time'] = file_parts[6]
        elif file['index'].startswith('noise'):
            file['type'] = 'noise'
            file['swath'] = file_parts[2].upper()
            file['polarization'] = file_parts[4].upper()
            file['start_time'] = file_parts[5]
            file['stop_time'] = file_parts[6]
        elif file['index'].startswith('rfi'):
            file['type'] = 'rfi'
            file['swath'] = file_parts[2].upper()
            file['polarization'] = file_parts[4].upper()
            file['start_time'] = file_parts[5]
            file['stop_time'] = file_parts[6]
        else:
            file['type'] = 'measurement'
            file['swath'] = file_parts[1].upper()
            file['polarization'] = file_parts[3].upper()
            file['start_time'] = file_parts[4]
            file['stop_time'] = file_parts[5]
        xml = (f'/dataObjectSection/dataObject[{i+1}]/byteStream/@size')
        file['size'] = int(doc.xpath(xml)[0])
        xml = (f'/dataObjectSection/dataObject[{i+1}]/byteStream/checksum')
        file['checksum'] = doc.xpath(xml)[0].text
        files.append(file)

    return files
