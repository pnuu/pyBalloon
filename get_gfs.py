#!/usr/bin/python
import sys
import requests
import datetime as dt

def get_gfs_data(datestr, utc_hour, area, verbose=False):
    """Download GFS data from NOMADS-NOAA for requested time and area.

    Required arguments:
        - datestr -- requested date as a string, eg. '20131803'
        - utc_hour -- requested UTC hour with two digits as a string
        - area -- 4-tuple of coordinates limiting the area to be retrieved.
    Coordinates are given as floats: (blat, tlat, llon, rlon), where
        - blat -- bottom latitude (southern) limit (float)
        - tlat -- top latitude (northern) limit (float)
        - llon -- left longitude (western) limit (float)
        - rlon -- right longitude (eastern) limit (float)

    Latitudes are north-positive, longitudes east-positive. Limits are
    inclusive. Data from closest available time will be downloaded for
    the latest available GFS main run, EPS main and 20 EPS members.
    """

    # Read requested date, time and  area

    yyyy, mm, dd = int(datestr[0:4]), int(datestr[4:6]), int(datestr[6:8])
    request_time = dt.datetime(yyyy, mm, dd, int(utc_hour))
    # 
    url = 'http://nomads.ncep.noaa.gov/cgi-bin/filter_gens.pl'
    req = requests.get(url)
    if req.error is not None:
        print "Could not connect! Error code: " % req.error
        sys.exit()

    text = req.content.split('</a>')
    available_days = []
    for t in text:
        if 'gefs' in t:
            available_days.append(t.split('gefs.')[-1])

    url_base = 'http://nomads.ncep.noaa.gov/cgi-bin/' \
        'filter_gens.pl?dir=%2Fgefs.'

    llon, rlon, llat, blat = area

    good_day = None
    good_init = None

    for day in available_days:
        yyyy2, mm2, dd2 = int(day[0:4]), int(day[4:6]), int(day[6:8])

        url = url_base+day
        req = requests.get(url)
        text = req.content.split('</a>')
        available_inits = []
        for t in text:
            if 'gefs' in t:
                available_inits.append(t.split('gefs.')[-1].split('>')[-1])

        for init in available_inits:
            # Calculate correct step for requested launch date/time
            hh2 = int(init)
            init_dt = dt.datetime(yyyy2, mm2, dd2, hh2)
            delta_t = request_time - init_dt

            if delta_t.seconds < 0 or delta_t.days < 0:
                continue
        
            # 00 06 12 18 ... for ensemble members
            ens_step = '%02d' % (6*int(delta_t.seconds/(6.*60*60)))
            # 00 03 06 09 12 ... for 0.5 and 1.0 deg main runs
            main_step = '%02d' % (3*int(delta_t.seconds/(3.*60*60)))
        
            for ens in range(1, 21):
                ens = '%02d' % ens

                ens_url = 'http://nomads.ncep.noaa.gov/cgi-bin/' \
                    'filter_gens.pl?file=gep'+ens+'.t'+init+ \
                    'z.pgrb2f'+ens_step+'&lev_1000_mb=on' \
                    '&lev_100_mb=on&lev_10_mb=on&' \
                    'lev_150_mb=on&lev_200_mb=on&lev_20_mb=on&' \
                    'lev_250_mb=on&lev_2_m_above_ground=on&' \
                    'lev_300_mb=on&lev_30_mb=on&lev_350_mb=on&' \
                    'lev_400_mb=on&lev_450_mb=on&lev_500_mb=on&' \
                    'lev_50_mb=on&lev_550_mb=on&lev_600_mb=on&' \
                    'lev_650_mb=on&lev_700_mb=on&lev_70_mb=on&' \
                    'lev_750_mb=on&lev_800_mb=on&lev_850_mb=on&' \
                    'lev_900_mb=on&lev_925_mb=on&lev_950_mb=on&' \
                    'lev_975_mb=on&lev_surface=on&' \
                    'var_HGT=on&var_TMP=on&var_UGRD=on&var_VGRD=on&' \
                    'subregion=&leftlon='+str(llon)+'&rightlon='+ \
                    str(rlon)+'&toplat='+str(tlat)+'&bottomlat='+str(blat)+ \
                    '&dir=%2Fgefs.'+day+'%2F'+init+'%2Fpgrb2'
            
                if verbose:
                    print ens_url

                ens_out = 'ens_' + ens + '.grib2'
                req = requests.get(ens_url)

                if req.status_code != 200:
                    print "Could not get ensemble data!"
                    break

                print "Saving ensemble member %d" % int(ens)
                fid = open(ens_out, 'wb')
                fid.write(req.content)
                fid.close()
            
                good_day = day
                good_init = init

                if ens == '20':
                    break

            if good_day is not None:
                break

        if good_day is not None:
            break

    day = good_day
    init = good_init
    step = main_step

    ens_main_url = 'http://nomads.ncep.noaa.gov/cgi-bin/filter_gfs.pl?file=' \
        'gfs.t'+init+'z.pgrbf'+step+'.grib2&lev_1000_mb=on&lev_100_mb=on' \
        '&lev_10_mb=on&lev_150_mb=on&lev_200_mb=on&lev_20_mb=on&lev_250_mb=on' \
        '&lev_2_m_above_ground=on&lev_300_mb=on&lev_30_mb=on&lev_350_mb=on' \
        '&lev_400_mb=on&lev_450_mb=on&lev_500_mb=on&lev_50_mb=on' \
        '&lev_550_mb=on&lev_600_mb=on&lev_650_mb=on&lev_700_mb=on' \
        '&lev_70_mb=on&lev_750_mb=on&lev_800_mb=on&lev_850_mb=on' \
        '&lev_900_mb=on&lev_925_mb=on&lev_950_mb=on&lev_975_mb=on' \
        '&lev_surface=on&var_HGT=on&var_TMP=on&var_UGRD=on&var_VGRD=on' \
        '&subregion=&leftlon='+str(llon)+'&rightlon='+str(rlon)+ \
        '&toplat='+str(tlat)+'&bottomlat='+str(blat)+'&dir=%2Fgfs.'+ \
        day+init

    if verbose:
        print ens_main_url

    req = requests.get(ens_main_url)
    if req.status_code != 200:
        print "Could not get ensemble main data!"
        sys.exit()
    print "Saving ensemble main"
    fid = open('ens_main.grib2', 'wb')
    fid.write(req.content)
    fid.close()

    main_url = 'http://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_hd.pl?' \
        'file=gfs.t'+init+'z.mastergrb2f'+step+'&lev_1000_mb=on' \
        '&lev_100_mb=on&lev_10_mb=on&lev_125_mb=on&lev_150_mb=on' \
        '&lev_175_mb=on&lev_1_mb=on&lev_20_mb=on&lev_225_mb=on' \
        '&lev_250_mb=on&lev_275_mb=on&lev_2_m_above_ground=on' \
        '&lev_2_mb=on&lev_300_mb=on&lev_30_mb=on&lev_325_mb=on' \
        '&lev_350_mb=on&lev_375_mb=on&lev_3_mb=on&lev_400_mb=on' \
        '&lev_425_mb=on&lev_450_mb=on&lev_475_mb=on&lev_500_mb=on' \
        '&lev_525_mb=on&lev_550_mb=on&lev_575_mb=on&lev_5_mb=on' \
        '&lev_600_mb=on&lev_625_mb=on&lev_650_mb=on&lev_675_mb=on' \
        '&lev_700_mb=on&lev_70_mb=on&lev_725_mb=on&lev_750_mb=on' \
        '&lev_775_mb=on&lev_7_mb=on&lev_800_mb=on&lev_825_mb=on' \
        '&lev_850_mb=on&lev_875_mb=on&lev_900_mb=on&lev_925_mb=on' \
        '&lev_950_mb=on&lev_975_mb=on&lev_surface=on&var_HGT=on' \
        '&var_TMP=on&var_UGRD=on&var_VGRD=on&subregion=&leftlon='+ \
        str(llon)+'&rightlon='+str(rlon)+'&toplat='+str(tlat)+ \
        '&bottomlat='+str(blat)+'&dir=%2Fgfs.'+day+init+'%2Fmaster'

    if verbose:
        print main_url

    print "Saving GFS main"
    req = requests.get(main_url)
    fid = open('gfs_main.grib2', 'wb')
    fid.write(req.content)
    fid.close()

    print '\n', 'Retrieved data:'
    print 'GFS main run', good_day, good_init+'Z +', main_step, 'h'
    print 'Ensemble main run', good_day, good_init+'Z +', main_step, 'h'
    print 'GFS main and EPS main valid time', datestr, \
        '%02dZ' % (int(good_init)+int(main_step))
    print "Ensemble members' run", good_day, good_init+'Z +', ens_step, 'h'
    print 'Ensemble valid time', \
        datestr, '%02dZ' % (int(good_init)+int(ens_step))


if __name__ == '__main__':
    """Command-line interface for downloading GFS data from
    NOMADS-NOAA for requested time and area.

    Usage:
    python get_gfs.py <yyyymmdd> <utc_hour> <blat,tlat,llon,rlon>

    where
        - yyyymmdd -- requested date
        - utc_hour -- requested UTC hour with two digits
        - blat -- bottom latitude (southern) limit
        - tlat -- top latitude (northern) limit
        - llon -- left longitude (western) limit
        - rlon -- right longitude (eastern) limit

    Latitudes are north-positive, longitudes east-positive. Limits are
    inclusive. Data from closest available time will be downloaded for
    the latest available GFS main run, EPS main and 20 EPS members.

    Example:
    python get_gfs.py 20130318 09 59.0,62.0,24.5,27.5
    """

    datestr = sys.argv[1]
    utc_hour = sys.argv[2]
    blat, tlat, llon, rlon = sys.argv[3].split(',')
    area = (float(blat), float(tlat), float(llon), float(rlon))

    get_gfs_data(datestr, utc_hour, area, verbose=True)






'''
rewrite-stuff, not yet completed

def get_available_days(run_type='main'):
    '''
    '''
    srch = 'gfs'
    url = 'http://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_hd.pl'

    if run_type == 'ens':
        srch = 'gefs'
        url = 'http://nomads.ncep.noaa.gov/cgi-bin/filter_gens.pl'

    while True:
        try:
            req = requests.get(url)
            if req.error is not None:
                print "Could not connect! Error code: " % req.error
                sys.exit()
            break
        except:
            print 'Request failed, retrying: ' + url

    text = req.content.split('</a>')
    available_days = []
    for t in text:
        if srch in t:
            available_days.append(t.split(srch+'.')[-1])

    available_dates = []
    for d in available_days:
        if run_type == 'main':
            available_dates.append(dt.datetime(int(d[:4]), 
                                               int(d[4:6]), 
                                               int(d[6:8]), 
                                               int(d[8:]), 0, 0))
        else:
            available_dates.append(dt.datetime(int(d[:4]), 
                                               int(d[4:6]), 
                                               int(d[6:]), 
                                               0, 0, 0)) # 00 UTC

    return available_dates


def get_ens_inits(date):
    
    base_url = 'http://nomads.ncep.noaa.gov/cgi-bin/' \
        'filter_gens.pl?dir=%2Fgefs.'
    srch = 'gefs'

    while True:
        try:
            url = base_url + '%4d%02d%02d' % (date.year, 
                                              date.month, 
                                              date.day)

            req = requests.get(url)

            if req.error is not None:
                print "Could not connect! Error code: " % req.error
                sys.exit()
            break
        except:
            print 'Request failed, retrying: ' + url

    text = req.content.split('</a>')
    available_inits = []
    for t in text:
        if srch in t:
            hh = int(t.split(srch+'.')[-1].split('>')[-1])
            available_inits.append(dt.datetime(date.year,
                                               date.month,
                                               date.day,
                                               hh))

    return available_inits


def download(date, step, area, run_type="main", ens_member=None):

    blat, tlat, llon, rlon = area
    day = "%4d%02d%02d" % (date.year, date.month, date.day)
    init = "%02d" % (date.hour)
    
    step = '%02d' % step

    url = 'http://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_hd.pl?' \
        'file=gfs.t'+init+'z.mastergrb2f'+step+'&lev_1000_mb=on' \
        '&lev_100_mb=on&lev_10_mb=on&lev_125_mb=on&lev_150_mb=on' \
        '&lev_175_mb=on&lev_1_mb=on&lev_20_mb=on&lev_225_mb=on' \
        '&lev_250_mb=on&lev_275_mb=on&lev_2_m_above_ground=on' \
        '&lev_2_mb=on&lev_300_mb=on&lev_30_mb=on&lev_325_mb=on' \
        '&lev_350_mb=on&lev_375_mb=on&lev_3_mb=on&lev_400_mb=on' \
        '&lev_425_mb=on&lev_450_mb=on&lev_475_mb=on&lev_500_mb=on' \
        '&lev_525_mb=on&lev_550_mb=on&lev_575_mb=on&lev_5_mb=on' \
        '&lev_600_mb=on&lev_625_mb=on&lev_650_mb=on&lev_675_mb=on' \
        '&lev_700_mb=on&lev_70_mb=on&lev_725_mb=on&lev_750_mb=on' \
        '&lev_775_mb=on&lev_7_mb=on&lev_800_mb=on&lev_825_mb=on' \
        '&lev_850_mb=on&lev_875_mb=on&lev_900_mb=on&lev_925_mb=on' \
        '&lev_950_mb=on&lev_975_mb=on&lev_surface=on&var_HGT=on' \
        '&var_TMP=on&var_UGRD=on&var_VGRD=on&subregion=&leftlon='+ \
        str(llon)+'&rightlon='+str(rlon)+'&toplat='+str(tlat)+ \
        '&bottomlat='+str(blat)+'&dir=%2Fgfs.'+day+init+'%2Fmaster'
    out = 'gfs_main.grib2'

    if run_type == 'ens':
        ens_member = '%02d' % ens_member
        url = 'http://nomads.ncep.noaa.gov/cgi-bin/' \
            'filter_gens.pl?file=gep'+ens_member+'.t'+init+ \
            'z.pgrb2f'+step+'&lev_1000_mb=on' \
            '&lev_100_mb=on&lev_10_mb=on&' \
            'lev_150_mb=on&lev_200_mb=on&lev_20_mb=on&' \
            'lev_250_mb=on&lev_2_m_above_ground=on&' \
            'lev_300_mb=on&lev_30_mb=on&lev_350_mb=on&' \
            'lev_400_mb=on&lev_450_mb=on&lev_500_mb=on&' \
            'lev_50_mb=on&lev_550_mb=on&lev_600_mb=on&' \
            'lev_650_mb=on&lev_700_mb=on&lev_70_mb=on&' \
            'lev_750_mb=on&lev_800_mb=on&lev_850_mb=on&' \
            'lev_900_mb=on&lev_925_mb=on&lev_950_mb=on&' \
            'lev_975_mb=on&lev_surface=on&' \
            'var_HGT=on&var_TMP=on&var_UGRD=on&var_VGRD=on&' \
            'subregion=&leftlon='+str(llon)+'&rightlon='+ \
            str(rlon)+'&toplat='+str(tlat)+'&bottomlat='+str(blat)+ \
            '&dir=%2Fgefs.'+day+'%2F'+init+'%2Fpgrb2'

        out = 'ens_' + ens_member + '.grib2'

    if run_type == 'ens_main':
        url = 'http://nomads.ncep.noaa.gov/cgi-bin/filter_gfs.pl?file=' \
            'gfs.t'+init+'z.pgrbf'+step+'.grib2&lev_1000_mb=on&lev_100_mb=on' \
            '&lev_10_mb=on&lev_150_mb=on&lev_200_mb=on&lev_20_mb=on&lev_250_mb=on' \
            '&lev_2_m_above_ground=on&lev_300_mb=on&lev_30_mb=on&lev_350_mb=on' \
            '&lev_400_mb=on&lev_450_mb=on&lev_500_mb=on&lev_50_mb=on' \
            '&lev_550_mb=on&lev_600_mb=on&lev_650_mb=on&lev_700_mb=on' \
            '&lev_70_mb=on&lev_750_mb=on&lev_800_mb=on&lev_850_mb=on' \
            '&lev_900_mb=on&lev_925_mb=on&lev_950_mb=on&lev_975_mb=on' \
            '&lev_surface=on&var_HGT=on&var_TMP=on&var_UGRD=on&var_VGRD=on' \
            '&subregion=&leftlon='+str(llon)+'&rightlon='+str(rlon)+ \
            '&toplat='+str(tlat)+'&bottomlat='+str(blat)+'&dir=%2Fgfs.'+ \
            day+init
        out = 'ens_main.grib2'

    while True:
        try:
            req = requests.get(url)
        except:
            print "Connection failed, trying again."

        try:
            if req.error is None:
                break
        except:
            pass

    print "Saving %s" % out
    fid = open(out, 'wb')
    fid.write(req.content)
    fid.close()
    print "File saved."


def find_closest_time_and_step(wanted_time, available_times, run_type='main'):
    '''
    '''

    if run_type == 'main':
        hour_step = 3
    else:
        hour_step = 6

    # Take the newest run as baseline
    closest = available_times[0]
    closest_diff = closest - wanted_time

    if closest_diff.total_seconds()/3600. >= hour_step/2:
        closest = available_times[1]
        closest_diff = closest - wanted_time
        for i in xrange(2, len(available_times)):
            td = available_times[i] - wanted_time
            if abs(td) <= abs(closest_diff):
                closest_diff = td
                closest = available_times[i]

    step = int(abs(round((closest_diff.total_seconds()/(3600*hour_step)))))
        
    return closest, step

'''
