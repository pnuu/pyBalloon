#!/usr/bin/python
"""Example main for pyBalloon"""

import numpy as np
import sys
import time
import pyb_io
import pyb_traj

#import matplotlib as mpl
#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap

if __name__ == "__main__":

    in_dir = sys.argv[1]

    time0 = time.time()

    # example balloon: 800 g -> http://www.novalynx.com/400-balloons.html
    balloon = {}
    balloon['altitude_step'] = 100.0 # meters
    balloon['equip_mass'] = 1.150 # kg
    balloon['balloon_mass'] = 0.8 # kg
    balloon['fill_radius'] = 0.925 #1.85/2 # meters
    balloon['radius_empty'] = 1.08/2 # meters
    balloon['burst_radius'] = 6.8/2 # meters
    balloon['thickness_empty'] = 0.2 * 10**-3 # mm -> meters
    balloon['Cd_balloon'] = 0.47
    balloon['Cd_parachute'] = 0.8
    balloon['parachute_areas'] = np.pi * np.array([0.5, 1.5])**2 # m^2
    balloon['parachute_change_altitude'] = 2000. # meters

    #num_random = 1000

    lat0 = 60.1 # degrees
    lon0 = 25.0 # degrees
    alt0 = 10. # meters

    model_data = pyb_io.read_gfs_set(in_dir, 
                                    (lat0+1.5, lon0-1.5, 
                                     lat0-0.5, lon0+1.5))#,
#                                     ens_main=None,
#                                     ens_member_pattern=None, alt0=alt0)

    print 'GFS data read, %.1f s elapsed' % (time.time() - time0)

    loc0 = (lat0, lon0, alt0)

    trajectories = []

    for data in model_data:
        trajectories.append(pyb_traj.calc_movements(data, loc0, balloon))

    print 'Trajectories calculated, %.1f s elapsed' % (time.time() - time0)

    # highest point in main-run trajectory
    idx, = np.where(trajectories[0]['alts'] == np.max(trajectories[0]['alts']))
    latx, _ = trajectories[0]['lats'][idx]
    lonx, _ = trajectories[0]['lons'][idx]
    altx, _ = trajectories[0]['alts'][idx]
    timex, _ = trajectories[0]['times'][idx]
    print latx, lonx, altx, '%.0f minutes' % (timex)
    other_info = [(latx, lonx, altx, 'Burst point', '%.0f minutes, %.0f meters' % (timex, altx))]

    kml_fname = '/tmp/pyballoon_trajectories.kml'
    pyb_io.save_kml(kml_fname, trajectories, other_info=other_info)
    
    print 'Program finished in %.1f s' % (time.time() - time0)

    """
    while True:
        # Read real balloon trajectory if available
        if live_data_fname is not None:
            live_data = read_live_data(live_data_fname)

        trajectories = []
        for data in model_data:
            trajectories.append(calc_movements(data, 
                                               loc0, 
                                               balloon, 
                                               live_data=live_data))

        save_kml(kml_fname, trajectories)
        time.sleep(60)
    """

    """
    print np.min(data['u_winds_std']), np.max(data['u_winds_std']), \
        np.min(data['v_winds_std']), np.max(data['v_winds_std'])
    """

    """
    np.savez('result.npz', all_lats=all_lats, all_lons=all_lons, 
             all_alts=all_alts, total_time=total_time, \
             distance_travelled=distance_travelled)
    """

    """
    for i in range(0, num_random+1):
        print "%3.1f %3.1f %2.3f %2.3f" % (
            total_time[-1], distance_travelled[i],
            all_lats[-1, i], all_lons[-1, i])
    """

    """
    t = total_time
    fig = plt.figure()
    ax = fig.gca()
    ax.set_aspect('equal')
    plt.plot(all_lons, all_lats, 'b')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.show()
    """
    '''
    fig = plt.figure()
    m = Basemap(projection='stere', lon_0=25, lat_0=60, lat_ts=60, 
                llcrnrlat=59.9, urcrnrlat=61.5, 
                llcrnrlon=24.5, urcrnrlon=26,
                rsphere=6371200.,resolution='l',area_thresh=10000)
    m.drawcoastlines()
    m.drawcountries()
    i = 0
    trajectories.reverse()
    for t in trajectories:
        x, y = m(t['lons'], t['lats'])
        if i == len(trajectories)-1:
            plt.plot(x, y, 'r')
        else:
            plt.plot(x, y, 'k.')
        i += 1
    #plt.ylabel('Latitude (deg)')
    #plt.xlabel('Longitude (deg)')
    #plt.legend(['U','V','Total'])
    plt.show()
    '''

