"""Functions used by pyBalloon to calculate balloon trajectories"""

import numpy as np
import pyb_aux

def calc_movements(data, loc0, balloon, live_data=None, alt_change_fit=10):
    """Calculate needed data and trajectory for a balloon. Data from
    live balloon flight can be used to replace past model values. In
    this case, the trajectory will be calculated starting from the
    last known real location.

    Required arguments:
        - data -- dictionary containing atmospheric data
        - loc0 -- initial location of the balloon
        - balloon -- dictionary containing balloon parameters
    
    Optional arguments:
        - live_data -- live data from balloon with 'lats', 'lons', 'alts'
        - alt_change_fit -- number of points to use for determining if the 
        - balloon is going up or down. Default: 10

    Return:
        - Dictionary containing calculated trajectory of the balloon
        containing 'lats', 'lons', 'alts', 'times' and 'distance'

    Todo:

    If live data are present:
    * Primary features re: live data
        - set existing trajectory from live data
        - calculate movement starting from latest known position:
    lats_live[-1], lons_live[-1], alts_live[-1]
    * Secondary features re: live data
        - if pressure and/or temperature are known, use (also) those 
        values for (initial? initial and descent?) lift calculations 
        (aaargh?)
        - calculate wind components from live data, use for descent
        - same for existing pressure/temperature data
        - or adjust model data with this data?
    """

    lat0, lon0, alt0 = loc0

    if live_data is not None:
        lats_live = live_data['lats']
        lons_live = live_data['lons']
        alts_live = live_data['altitudes']
        pres_live = live_data['pressures']
        temps_live = live_data['temperatures']

        # Check balloon direction
        # Default to "going up"
        if len(alts_live) > alt_change_fit:
            fit = np.polyfit(np.arange(0, alt_change_fit), 
                             alts_live[-alt_change_fit-1:-1], 1)
            if fit[0] < 0:
                # Balloon is going down
                pass

    data['air_densities'] = pyb_aux.air_density(data)
    data['balloon_radii'], gas_mass = \
        pyb_aux.mooney_rivlin(data, 
                              balloon['radius_empty'], 
                              balloon['fill_radius'], 
                              balloon['thickness_empty'])
    data['balloon_volumes'] = pyb_aux.balloon_volume(data)

    total_mass = balloon['equip_mass'] + \
        balloon['balloon_mass'] + gas_mass # kg

    data['lifts'] = pyb_aux.lift(data, total_mass)

    data['ascent_speeds'] = pyb_aux.ascent_speed(data, total_mass, 
                                                balloon['Cd_balloon'])

    data = pyb_aux.data_interpolation(data, 
                                      alt0, 
                                      balloon['altitude_step'], 
                                      mode='spline')

    data['descent_speeds'] = \
        pyb_aux.descent_speed(data, 
                             balloon['equip_mass'],
                             balloon['Cd_parachute'],
                             balloon['parachute_areas'], 
                             balloon['parachute_change_altitude'])

    max_altitudes, max_alt_idxs = \
        pyb_aux.burst_altitude(data, balloon['burst_radius'])

    # fix ...
    #zero_lift_alts, zero_lift_idxs = pyb_aux.neutral_buoyancy_level(data)
    # ... this to return arrays, if needed

    delta_t = balloon['altitude_step'] / data['ascent_speeds']
    # Mark levels above neutral buoyancy as invalid
    # not really needed
    '''
    idxs = np.where(delta_t < 0)
    #if(len(idxs[0]) > 0):
        idx = idxs[0][0]
        delta_t[idx:] = np.nan
    '''
    data['ascent_time_steps'] = delta_t
    data['cumulative_ascent_times'] = np.cumsum(delta_t)/60.

    delta_t = -1*balloon['altitude_step'] / data['descent_speeds']
    data['descent_time_steps'] = delta_t
    data['cumulative_descent_times'] = np.cumsum(delta_t)/60

    alts = data['altitudes']

    dxs_up = data['u_winds'] * data['ascent_time_steps']
    dys_up = data['v_winds'] * data['ascent_time_steps']

    dxs_down = data['u_winds'] * data['descent_time_steps']
    dys_down = data['v_winds'] * data['descent_time_steps']

    data_lats = np.radians(data['lats'])
    data_lons = np.radians(data['lons'])
    lat_rad = [np.radians(lat0)]
    lon_rad = [np.radians(lon0)]
    all_alts = [alt0]
    total_time = [0]
    distance_travelled = 0

    i = 0
    # Calculate the movement during ascent
    while True:
        # Find the closest grid point
        diff = np.sqrt((data_lats-lat_rad[-1])**2 +
                       (data_lons-lon_rad[-1])**2)
        grid_i, = np.where(diff == diff.min())
        grid_i = grid_i[0]
    
        lat, lon, dist = movement2ll(lat_rad[-1], 
                                     lon_rad[-1], 
                                     alts[i], 
                                     dxs_up[grid_i, i], 
                                     dys_up[grid_i, i])
        lat_rad.append(lat)
        lon_rad.append(lon)
        total_time.append(data['ascent_time_steps'][grid_i, i])
        all_alts.append(alts[i])
        distance_travelled += dist

        if all_alts[-1] >= max_altitudes[grid_i]:
            break

        i += 1

    # Calculate the movement during descent
    while i >= 0:
        # Find the closest grid point
        diff = np.sqrt((data_lats-lat_rad[-1])**2 +
                       (data_lons-lon_rad[-1])**2)
        grid_i, = np.where(diff == diff.min())
        grid_i = grid_i[0]

        lat, lon, dist = movement2ll(lat_rad[-1], \
                                         lon_rad[-1], \
                                         alts[i], \
                                         dxs_down[grid_i, i], \
                                         dys_down[grid_i, i])
        lat_rad.append(lat)
        lon_rad.append(lon)
        total_time.append(data['descent_time_steps'][grid_i, i])
        all_alts.append(alts[i])
        distance_travelled += dist

        i -= 1

    # Convert the result array lists to Numpy 2D-arrays
    output = {}
    output['lats'] = np.degrees(np.array(lat_rad)) # to decimal degrees
    output['lons'] = np.degrees(np.array(lon_rad)) # to decimal degrees
    output['alts'] = np.array(all_alts)
    output['times'] = np.cumsum(np.array(total_time))/60 # to minutes
    output['distance'] = distance_travelled

    print("Maximum altitude", np.max(all_alts), \
        'm, distance travelled %.1f:' % distance_travelled, 'km')
    print('Landing location', \
        '%.3f, %.3f' % (output['lats'][-1], output['lons'][-1]), \
        'after %d min of flight\n' % (int(output['times'][-1])))

    return output


def movement2ll(lat_rad, lon_rad, alt, dx, dy):
    """Calculate new lat/lon coordinates from Cartesian dx (east
    positive) and dy (north positive) displacements. Calculation is
    done using Haversine formula with radius derived from WGS84
    reference ellipsoid with altitude added to it.

    Required arguments:
        - lat_rad -- Current latitude in radians
        - lon_rad -- Current longitude in radians
        - alt -- Current altitude, in meters
        - dx -- East - west movement in meters (east positive)
        - dy -- North - south movement in meters (north positive)

    Return:
        - New coordinates: latitude [radians], longitude [radians],
        and distance traveled [km]
    """
    
    Re = pyb_aux.earth_radius(lat_rad)
    radius = Re + alt/1000. # Altitude is given in meters, convert to km
    
    # Distance
    dist = np.sqrt(dx*dx + dy*dy)/1000. # Convert to km
    # Direction
    theta = np.arctan2(dx, dy)

    cos_dr = np.cos(dist/radius)
    sin_dr = np.sin(dist/radius)
    sin_lat = np.sin(lat_rad)
    cos_lat = np.cos(lat_rad)
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)

    lat2 = np.arcsin(sin_lat * cos_dr + cos_lat * sin_dr * cos_theta)

    lon2 = lon_rad + np.arctan2(sin_theta * sin_dr * cos_lat, 
                                cos_dr - sin_lat * np.sin(lat2))

    return lat2, lon2, dist

