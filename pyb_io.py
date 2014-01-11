"""Input and output functions used by pyBalloon"""
import pygrib as pg
import numpy as np
import os
from glob import glob
import csv
import pyb_aux
import requests

def read_gfs_file(fname, area=None, alt0=0, t_0=None, extra_data=None):
    """Collect relevant information from GFS model.

    Required arguments:
        - fname -- File to read

    Optional arguments:
        - area -- tuple of NW and SE limits of the area to be read,
        eg.  (62., 22., 59., 24.). Default: None (read all data)
        - alt0 -- Starting altitude above sea level, default: 0.0 m
        - t_0 -- Temperature at ground level, default: None
        - extra_data -- Add highest levels from extra_data to GFS
        ensembles which are limited to 10 hPa.

    Return the following data for the closest pixel location as Numpy
    arrays in a dictionary:
        - u_winds -- U component of wind [m/s] - E-W component,
        positive *towards* East
        - v_winds -- V component of wind [m/s] - N-S component,
        positive *towards* North
        - temperatures -- Temperature [K]
        - altitudes -- Geopotential height [m]
    """
    
    if area is not None:
        tlat, llat, blat, rlat = area
    else:
        tlat, llat, blat, rlat = 90., -180., -90., 180.

    grib = pg.open(fname)
    grib.seek(0)
    u_msgs = grib.select(name='U component of wind')
    v_msgs = grib.select(name='V component of wind')
    g_msgs = grib.select(name='Geopotential Height')
    t_msgs = grib.select(name='Temperature')

    lats, lons = u_msgs[0].latlons()

    # Find closest pixel location
    locs = pyb_aux.all_and([lats <= tlat, lats >= blat, 
                           lons <= rlat, lons >= llat])
    row_idx, col_idx = np.where(locs)
    lats = lats[row_idx, col_idx]
    lons = lons[row_idx, col_idx]

    # Collect U component of wind data
    u_wind = {}
    for msg in u_msgs:
        if msg.typeOfLevel == 'isobaricInhPa':
            u_wind[msg.level] = msg.values[row_idx, col_idx]

    # Collect V component of wind data
    v_wind = {}
    for msg in v_msgs:
        if msg.typeOfLevel == 'isobaricInhPa':
            v_wind[msg.level] = msg.values[row_idx, col_idx]

    # Collect temperatures
    temperature = {}
    for msg in t_msgs:
        if msg.typeOfLevel == 'isobaricInhPa':
            temperature[msg.level] = msg.values[row_idx, col_idx]
        # Add msg.typeOfLevel == 'surface', save to another variable for
        # later use
        if msg.typeOfLevel == 'surface':
            t_surface = msg.values[row_idx, col_idx]

    # Collect Geopotential heights
    altitude = {}
    for msg in g_msgs:
        if msg.typeOfLevel == 'isobaricInhPa':
            altitude[msg.level] = msg.values[row_idx, col_idx]

    # Collect data to correct altitude order. Set "surface" values
    # before real data.
    u_winds = [np.zeros(lats.shape)]
    v_winds = [np.zeros(lats.shape)]
    # Use given surface temperature if available, otherwise use the
    # model value
    if t_0 is None:
        temperatures = [t_surface]
    else:
        temperatures = [t_0*np.ones(lats.shape)]
    altitudes = [alt0*np.ones(lats.shape)]
    pressures = u_wind.keys()
    pressures.append(max(pressures))

    # Put pressures in altitude order and use them as keys
    pressures.sort()
    pressures.reverse()
    i = 0
    for key in pressures:
        if i == 0:
            i = 1
        else:
            uwnd, vwnd, temp, alt = [], [], [], []
            uwnd.append(u_wind[key])
            vwnd.append(v_wind[key])
            temp.append(temperature[key])
            alt.append(altitude[key])

            # Add extra data to complement the currently read data,
            # ie. 1, 2, 3, 5 and 7 hPa levels from GFS main run to
            # ensembles. Data are expected to be in the same format as
            # returned by this function.
            j = 0
            if extra_data is not None:
                for level in extra_data['pressures'][:, 0]/100:
                    if level < min(pressures):
                        for idx in range(0, len(extra_data['u_winds'][0, :])):
                            uwnd.append(extra_data['u_winds'][j, idx])
                            vwnd.append(extra_data['v_winds'][j, idx])
                            temp.append( \
                                extra_data['temperatures'][j, idx])
                            alt.append(extra_data['altitudes'][j, idx])
                    j += 1

            u_winds.append(np.hstack(uwnd))
            v_winds.append(np.hstack(vwnd))
            temperatures.append(np.hstack(temp))
            altitudes.append(np.hstack(alt))

    # Convert data in lists to Numpy arrays and add them to a
    # dictionary that is returned
    data = {}
    data['lats'] = np.array(lats)
    data['lons'] = np.array(lons)
    data['u_winds'] = np.array(u_winds)
    data['v_winds'] = np.array(v_winds)
    data['temperatures'] = np.array(temperatures)
    data['altitudes'] = np.array(altitudes)
    all_pressures = []
    for dat in data['lats']:
        all_pressures.append(100*np.array(pressures)) # Convert hPa to Pa
    data['pressures'] = np.array(all_pressures).transpose()

    return data


def read_gfs_set(directory, area=None, alt0=0, main='gfs_main.grib2',
                 ens_main='ens_main.grib2',
                 ens_member_pattern='ens_??.grib2', 
                 use_extra=False):
    """Read a set of GFS data consisting of 0.5 degree main run, 1
    degree ensemble main and 1 degree ensemble members. The 1 degree
    data are extended to 7 - 1 hPa levels with data from the main run.

    Required arguments:
        - directory -- Directory containing the data

    Optional arguments:
        - area -- tuple of NW and SE limits of the area to be read,
        eg.  (62., 22., 59., 24.). Default: None (read all data)
        - alt0 -- Starting altitude of the balloon, meters from the WGS84
        reference ellipsoid. Default: 0.0.
        - main -- Filename of the GFS main run. Default: gfs_main.grib2.
        - ens_main -- Filename of the GFS ensemble main run. Default:
        ens_main.grib2.
        - ens_member_pattern -- Filename pattern of the ensemble runs.
        Default: ens_??.grib2

    Return:
        - List of dictionaries containing u_winds, v_winds, temperatures,
        pressures and altitudes.
    """

    all_data = []

    fname = os.path.join(directory, main)
    print "Reading GFS operational run from", fname
    main_run_data = read_gfs_file(fname, area=area, alt0=alt0)
    all_data.append(main_run_data)

    if ens_main is not None:
        fname = os.path.join(directory, ens_main)
        print "Reading GFS ensemble main run from", fname
        if use_extra:
            ens_main_data = read_gfs_file(fname, 
                                          area=area,
                                          alt0=alt0, 
                                          extra_data=main_run_data)
        else:
            ens_main_data = read_gfs_file(fname, 
                                          area=area,
                                          alt0=alt0, 
                                          extra_data=None)
        all_data.append(ens_main_data)

    if ens_member_pattern is not None:
        ens_files = glob(os.path.join(directory, ens_member_pattern))
        ens_files.sort()
    
        for fname in ens_files:
            print "Reading GFS ensemble member from", fname

            if use_extra:
                all_data.append(read_gfs_file(fname, 
                                              area=area, 
                                              alt0=alt0, 
                                              extra_data=main_run_data))
            else:
                all_data.append(read_gfs_file(fname, 
                                              area=area, 
                                              alt0=alt0, 
                                              extra_data=None))

    return all_data


def get_sounding(station_id, date, utc_hour):
    """Get sounding data using station ID, date and time.

    Data are retrived from Universio of Wyoming. Valid times depend on
    the station, normal sounding times are 00 and 12 UTC, or 06 and
    18 UTC.

    Required arguments:
        - station_id -- ID of the station: number or air-port code
        - date -- Date as 3-tuple (yyyy, mm, dd), eg. (2013, 3, 18)
        - utc_hour -- UTC hour of the requested sounding eg. 6

    Returns:
        - Dictionary containing station coordinates, 'u_winds',
          'v_wind', 'temperatures', 'altitudes' and 'pressures' as
          Numpy arrays
    """

    knots2ms = 0.514

    year, month, day = date[0], date[1], date[2]

    url = 'http://weather.uwyo.edu/cgi-bin/sounding?' + \
        'TYPE=TEXT%3ALIST&YEAR=' + str(year) + '&MONTH=' + \
        '%02d' % month + '&FROM=' + '%02d%02d' % (day, utc_hour) + \
        '&TO=' + '%02d%02d' % (day, utc_hour) + '&STNM=' + str(station_id)

    req = requests.get(url)
    text  = req.text
    station_lat = float(text.split('Station latitude: ')[-1].split('\n')[0])
    station_lon = float(text.split('Station longitude: ')[-1].split('\n')[0])
    data = text.split('<PRE>')[1].split('</PRE>')[0].split('\n')
    data = data[5:] # only the numerical rows
    
    pressures = []
    altitudes = []
    temperatures = []
    u_winds = []
    v_winds = []

    for row in data:
        nums = row.split()
        if len(nums) == 11:
            pressures.append(float(nums[0]))
            altitudes.append(float(nums[1]))
            temperatures.append(273.15+float(nums[2]))
        
            wdir = np.radians(float(nums[6]))
            wspd = float(nums[7])*knots2ms
        
            # Towards North and East are positive
            u_winds.append(-1*wspd*np.sin(wdir))
            v_winds.append(-1*wspd*np.cos(wdir))

    data = {}
    data['lats'] = np.array(station_lat)
    data['lons'] = np.array(station_lon)
    data['u_winds'] = np.array(u_winds)
    data['v_winds'] = np.array(v_winds)
    data['temperatures'] = np.array(temperatures)
    data['pressures'] = 100*np.array(pressures)
    data['altitudes'] = np.array(altitudes)

    return data


def read_live_data(fname, delimiter=',', temp_conv=(1, 0),
                   pressure_conv=(1, 0)):
    """Read data from live feed (via a file) from a balloon (or
    simulator). The data are assumed to be in CSV format with the
    following fields:
        - latitude, longitude, altitude
    or
        - latitude, longitude, altitude, pressure
    or
        - latitude, longitude, altitude, pressure, temperature

    Location data are compulsory (first 3 fields). Pressure and
    temperature fields can be empty, missing or in use.

    Required arguments:
        - fname -- Filename of the live data.

    Optional arguments:
        - delimiter -- CSV field delimiter. Default: ','
        - temp_conv -- Conversion factors (multiplier and offset)
        between temperature data and temperature in Kelvins. Default
        (Kelvin -> Kelvin): (1, 0).
        - pressure_conv -- Conversion factors (multiplier and offset)
        between
        - pressure data and pressure in Pascals. Default: (Pa -> Pa):
        (1, 0)

    Return:
        - None or dictionary of Numpy arrays containing the data, if
        any available.
    """

    try:
        fid = open(fname)
    except:
        return None

    reader = csv.reader(fid, delimiter=delimiter)

    lats = []
    lons = []
    altitudes = []
    pressures = []
    temperatures = []
    num = 0
    for row in reader:
        if len(row) == 3:
            lat, lon, alt = row
            pres, temp = None, None
        else:
            if len(row) == 4:
                lat, lon, alt, pres = row
                temp = None
            else:
                lat, lon, alt = row[0], row[1], row[2]
                pres, temp = row[3], row[4]
        
        lats.append(float(lat))
        lons.append(float(lon))
        altitudes.append(float(alt))
        pressures.append(float(pres))
        temperatures.append(float(temp))
        num += 1

    if num == 0:
        return None

    live_data = {}
    live_data['lats'] = np.array(lats)
    live_data['lats'] = np.array(lons)
    live_data['altitudes'] = np.array(altitudes)
    live_data['pressures'] = np.array(pressures) * pressure_conv[0] + \
        pressure_conv[1]
    live_data['temperatures'] = np.array(temperatures) * temp_conv[0] + \
        temp_conv[1]

    return live_data


def save_kml(fname, data, model_start_idx=0, 
             eps_mode='end', other_info=None):
    """Save given trajectories as KML. The first trajectory is assumed
    to be from GFS main run, the second from ensemble main and the
    rest from other ensemble members.

    Required arguments:
        - fname -- File to save the KML data to
        - data -- List of dictionaries containing 'lats', 'lons',
        'altitudes' and 'times.

    Optional arguments
        - model_start_idx -- Vector index to show where model
        trajectory starts. Trajectory before this is from live
        data. Default: 0
        - eps_mode -- How to present ensemble predictions of
        trajectories. Possible modes are 'end' and 'full'. Default:
        'end'
        - other_info -- Additional information to save into the KML
        path: eg.  other_info = [(lat, lon, alt, 'marker name', 
                                  'marker description')]

    Todo:
        - modelled current location of the balloon
        - latest live data location
        - important times / altitudes / ...
        - maximum modelled altitude / maximum gained altitude (balloon
        burst at T+xx:xx:xx, xx.x km)
        - different colors for model / live data parts of the trajectory?
            - change color of the trajectory at model_start_idx
            - or different paths?
    """

    kml_str = '<?xml version="1.0" encoding="UTF-8"?>\n'
    kml_str += '<kml xmlns="http://www.opengis.net/kml/2.2">\n'
    kml_str += '<Document>\n'
    kml_str += '<name>pyballoon trajectory</name>\n'
    kml_str += '<description>pyballoon trajectories</description>\n'

    kml_str += '<Style id="yellowLineGreenPoly">\n'
    kml_str += '<LineStyle>\n'
    kml_str += '<color>2000ffff</color>\n'
    kml_str += '<width>4</width>\n'
    kml_str += '</LineStyle>\n'
    kml_str += '<PolyStyle>\n'
#    kml_str += '<color>4000ff00</color>\n'
    kml_str += '<color>0000ff00</color>\n'
    kml_str += '</PolyStyle>\n'
    kml_str += '</Style>\n'

    kml_str += '<Style id="redLineGreenPoly">\n'
    kml_str += '<LineStyle>\n'
    kml_str += '<color>7f0000ff</color>\n'
    kml_str += '<width>4</width>\n'
    kml_str += '</LineStyle>\n'
    kml_str += '<PolyStyle>\n'
#    kml_str += '<color>4000ff00</color>\n'
    kml_str += '<color>0000ff00</color>\n'
    kml_str += '</PolyStyle>\n'
    kml_str += '</Style>\n'

    num = 0
    for dat in data:
        if num == 0 or eps_mode == 'full':
            kml_str += '<Placemark>\n'
            if num == 0:
                kml_str += '<name>GFS main</name>\n'
            else:
                if num == 1:
                    kml_str += '<name>GFS Ensemble main</name>\n'
                else:
                    kml_str += '<name>GFS Ensemble member %d</name>\n' % \
                        (num-1)
            if num == 0:
                kml_str += '<description>pyballoon trajectory</description>\n'
                kml_str += '<styleUrl>#redLineGreenPoly</styleUrl>\n'
            else:
                if eps_mode == 'full':
                    kml_str += '<description>pyballoon trajectory' \
                        '</description>\n'
                    kml_str += '<styleUrl>#yellowLineGreenPoly</styleUrl>\n'

            kml_str += '<LineString>\n'
            kml_str += '<extrude>1</extrude>\n'
            kml_str += '<tessellate>1</tessellate>\n'
            kml_str += '<altitudeMode>absolute</altitudeMode>\n'
            kml_str += '<coordinates>\n'

            t_prev = -2.0
            for i in range(0, len(dat['lats'])):
                if dat['times'][i] >= t_prev + 0.5:
                    t_prev = dat['times'][i]
                    kml_str += '%f,%f,%f\n' % \
                        (dat['lons'][i], dat['lats'][i], dat['alts'][i])

            kml_str += '</coordinates>\n'
            kml_str += '</LineString>\n'
            kml_str += '</Placemark>\n'
        num += 1

    # Add placemarks for the trajectory end-points
    num = 0
    for dat in data:
        kml_str += '<Placemark>\n'
        if num == 0:
            kml_str += '<name>GFS main run</name>\n'
            kml_str += '<description>End-point based on GFS ' \
                'main run</description>\n'
        else:
            kml_str += '<name># %d</name>\n' % (num-1)
            kml_str += '<description>End-point based on GFS ' \
                    'ensemble member %d</description>\n' % (num-1)

        kml_str += '<Point>\n'
        kml_str += '<coordinates>%f,%f</coordinates>\n' % \
            (dat['lons'][-1], dat['lats'][-1])
        kml_str += '</Point>\n'
        kml_str += '</Placemark>\n'
        num += 1

    # Add "other_info" places
    for dat in other_info:
        kml_str += '<Placemark>\n'
        kml_str += '<name>'+dat[3]+'</name>\n'
        kml_str += '<description>'+dat[4]+'</description>\n'
        kml_str += '<Point>\n'
        kml_str += '<coordinates>%f,%f,%f</coordinates>\n' % \
            (dat[1], dat[0], dat[2])
        kml_str += '</Point>\n'
        kml_str += '</Placemark>\n'

    kml_str += '</Document>\n'
    kml_str += '</kml>\n'

    fid = open(fname, 'w')
    fid.write(kml_str)
    fid.close()

