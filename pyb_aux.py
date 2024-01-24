"""Auxiliary functions used in pyBalloon"""

import numpy as np
from scipy import interpolate

g_0 = 9.80665 # Earth gravitational acceleration at surface
R_e = 6371009 # mean Earth radius in meters
R = 8.3144621 # Ideal gas constant
M_air = 0.0289644 # molar mass of air [kg/mol], altitude dependence
                  # not used here
M_helium = 4.002602
Cd_sphere = 0.47 # Drag coefficient for a sphere


def all_and(data):
    """Logical and for a list of arrays.
    
    Required arguments:
        - data -- list of Numpy boolean arrays

    Return:
        - result -- Logical and of all given arrays
    """
    result = data.pop()
    for d in data:
        result = np.logical_and(result, d)

    return result


def earth_radius(lat_rad):
    """Calculate Earth WGS84 based radius on a given latitude.

    http://en.wikipedia.org/wiki/Earth_radius#Radius_at_a_given_geodetic_latitude
    """

    # WGS84 reference ellipsoid
    a = 6378.137     # Earth equatorial radius, km
    b = 6356.7523142 # Earth polar radius, km
    
    cos_lat = np.cos(lat_rad)
    sin_lat = np.sin(lat_rad)

    r = np.sqrt(((a*a*cos_lat)**2 + (b*b*sin_lat)**2) / 
                ((a*cos_lat)**2 + (b*sin_lat)**2))

    return r


def air_density(data):
    """Calculate air density from pressure and temperature.
    
    Required argument:
        - data -- Dictionary returned by read_gfs_data()

    Return:
        - Air densities (as Numpy array) for each altitude step in data.
    """

    p = data['pressures']
    T = data['temperatures']

    if p.shape != T.shape:
        x, y = p.shape
        rho = []
        for i in range(0, y): 
            rho.append(np.array((p[:, i] * M_air)/(R * T)))
        rho = np.array(rho).transpose()
    else:
        rho = (p * M_air)/(R * T)

    return rho


def data_interpolation(data, alt0, step, mode='spline'):
    """Interpolate (and extrapolate in the low end, if mode='spline'
    is used) vertical data from alt0 to maximum level present in the
    data.

    Required arguments:
        - data -- Dictionary returned by read_gfs_data()
        - alt0 -- Starting altitude of the interpolation in meters from
        WGS84 reference ellipsoid
        - step -- Interpolation altitude step

    Optional arguments:
        - mode -- Interpolation method. Supported methods are 'spline' and
        'linear'. Default: 'spline'

    Return:
        - Interpolated version of the data.
    """

    altitudes = data['altitudes']

    new_data = {}

    new_data['altitudes'] = np.arange(alt0, altitudes.max(), step)
    new_data['lats'] = data['lats']
    new_data['lons'] = data['lons']
    for key in list(data.keys()):
        if key != 'altitudes' and key != 'lats' and key != 'lons':
            new_arr = []
            d = data[key]

            try:
                x, y = d.shape
            except ValueError:
                x = 0
                y, = d.shape

            if mode == 'spline':
                if x > 0:
                    for i in range(0, y):
                        ok_idxs = altitudes[:, i] >= alt0
                        tck = interpolate.splrep(altitudes[ok_idxs, i], 
                                                 d[ok_idxs, i])
                        new_arr.append(np.array(interpolate.splev( \
                                    new_data['altitudes'], tck)))
                else:
                    tck = interpolate.splrep(altitudes, d)
                    new_arr.append(np.array(interpolate.splev( \
                                new_data['altitudes'], tck)))

            else: # use linear interpolation
                # There's something wrong here:
                for i in range(0, y):
                    for i in range(0, len(d)):
                        tck = interpolate.interp1d(altitudes[:, i], d[:, i])
                        new_arr.append(tck(new_data['altitudes']))

            new_data[key] = np.array(new_arr)

    return new_data


def lift(data, mass):
    """Calculate effective lift (force, in Newtons) caused by the
    balloon.

    Required arguments:
        - data -- Dictionary containing 'altitudes' (meters), balloon
        'volumes' (m^3) and 'air_densities' (kg/m^3)
        - mass -- Mass of the whole balloon

    Return:
        - Resultant lift force.
    """

    h = data['altitudes']
    V_b = data['balloon_volumes']
    rho_air = data['air_densities']
    g = g_0 * (R_e / (R_e + h))**2 # Gravitational acceleration at height h

    F_lift = g * (V_b*rho_air - mass)

    return F_lift


def balloon_volume(data):
    """Calculate volume of a sphere.

    Required argument:
        - data -- Dictionary containing 'balloon_radii' (in meters)

    Return:
        - Volume of the balloon at each level in data
    """

    r = data['balloon_radii']
    V = 4/3. * np.pi * r**3

    return V


def balloon_volume_ideal_gas(data, gas_mass, gas_molar_mass=M_helium):
    """Calculate gas (~balloon) volume based on ideal gas law: pV =
    nRT.

    Required arguments:
        - data -- Dictionary returned by read_gfs_data()
        - gas_mass -- Mass of the gas in question

    Optional arguments:
        - gas_molar_mass -- Gas molar mask (g/mol). Default: 4.002602 (Helium)

    Return:
        - Gas volume at each level in data
    """

    m = gas_mass # 
    M = gas_molar_mass/1000. # default to Helium, convert to kg/mol

    T = data['temperatures']
    p = data['pressures'] # pressure in Pascals

    # gas volume without the balloon
    V = m*R*T/(M*p) # pV = nRT, n = m/M

    return V


def burst_altitude(data, burst_radius):
    """Find the altitude where balloon radius gets greater than the
    given burst radius.

    Required arguments:
        - data -- Dictionary containing 'altitudes' and corresponding
        'balloon_radii'
        - burst_radius -- Balloon burst radius

    Return:
        - Altitude of burst and corresponding array index
    """

    radii = data['balloon_radii']

    alt_out = []
    i_out = []
    for radius in radii:
        diff = np.abs(radius-burst_radius)
        idx = np.where(diff == diff.min())    
        i_out.append(idx[0][0])
        alt_out.append(data['altitudes'][idx])

    return np.array(alt_out), np.array(i_out)


def neutral_buoyancy_level(data):
    """Find the level of neutral buoyancy (or more precise, the level
    where effective lift is closest to zero).

    Required arguments:
        - data -- Dictionary containing 'altitudes' and corresponding 'lifts'

    Return:
        - Altitude of neutral buoyancy and corresponding array index
    """

    alt_out = []
    i_out = []
    idx = 0
    for lft in data['lifts']:
        lft = np.abs(lft)
        idx2 = np.where(lft == lft.min())
        i_out.append(idx2[0][0])
        alt_out.append(data['altitudes'][idx, idx2])
        idx += 1

    return np.array(alt_out), np.array(i_out)


def ascent_speed(data, mass, Cd=Cd_sphere):
    """Calculate the rate of ascent (in m/s) for the inflated balloon
    at given levels.

    Required arguments:
        - data -- Dictionary with 'altitudes' and corresponding
        'air_densities', 'balloon_radii' and 'balloon_volumes'
        - mass -- Full balloon mass (in kg)

    Optional arguments:
        - Cd -- Coefficient of drag. Default: 0.47 (sphere)

    Return:
        - Ascent speed (m/s) at every level in input data.
    """

    m = mass
    rho = data['air_densities']
    A = np.pi*data['balloon_radii']**2
    V = data['balloon_volumes']
    h = data['altitudes']

    g = g_0 * (R_e / (R_e + h))**2

    Fb = V*rho*g # byoyance
    Fg = m*g     # gravity
    F = Fb-Fg

    # Set negative buyoyancies to zero (we won't get there)
    idxs = np.where(F <= 0)
    if len(idxs[0]) > 0:
        idx = idxs[0][0]
        F[idx:] = 1e-30

    v = np.sqrt(2*F/(rho*Cd*A))

    return v


def descent_speed(data, mass, Cd, areas, change_alt=None):
    """Calculate the rate of descent for deflated (burst) balloon with
    1 or 2 different sized parachutes with given areas, change
    altitude and drag-coefficient.

    Required arguments:
        - data -- Dictionary with 'altitudes', and corresponding 'air_densities'
        - mass -- Mass of the payload + assumed remnants of the balloon
        - Cd -- Coefficients of drag for one or two parachutes in a tuple
        - areas -- Effective areas (in m^2) of one or two parachutes in a tuple

    Optional arguments:
        - change_alt -- Altitude where first parachute is changed to
        the second one. If None, only one parachute is used. Default:
        None

    Return:
        - Rate of descent (m/s) for each level in input data
    """
    
    m = mass
    h = data['altitudes']
    g = g_0 * (R_e / (R_e + h))**2 # Gravitational acceleration at height h

    if change_alt is not None:
        idxs = h < change_alt
        speeds = []
        for rho in data['air_densities']:
            v = np.sqrt(2*m*g/(rho*Cd*areas[0]))
            v[idxs] = np.sqrt(2*m*g[idxs]/(rho[idxs]*Cd*areas[1]))
            speeds.append(v)
    else:
        for rho in data['air_densities']:
            v = np.sqrt(2*m*g/(rho*Cd*areas))
            speeds.append(v)

    return -1*np.array(speeds)
    

def mooney_rivlin(data, radius_empty, radius_filled, 
                  thickness_empty, gas_molar_mass=M_helium):
    """Calculate balloon radii for given
    pressures/temperatures/initial conditions using inversion of
    Mooney-Rivlin equation.

    See description of the equation at: 
    http://www.zmatt.net/weather-balloon-physics/

    Required arguments:
        - data -- Dictionary containing 'pressures' and 'temperatures'
        - radius_empty -- Radius of the empty balloon
        - radius_filled -- Radius when filled at ground level
        - thickness_empty -- Balloon rubber initial thickness
    
    Optional arguments:
        - gas_molar_mass -- Molar mass (g/mol) of the gas used to fill
        the balloon. Default: 4.002602 (Helium)

    Return:
        - Radius of the balloon at each level in input data, and the
        mass of the gas.
    """

    r0 = radius_empty # in meters
    r1 = radius_filled # in meters
    t0 = thickness_empty # in meters
    M = gas_molar_mass / 1000. # convert to kg/mol
    p = data['pressures']
    p0 = p[0]
    T = data['temperatures']
    T0 = T[0]

    mu = 300000. # balloon shear modulus in Pascals
    alfa = 10./11.

    # Amount of gas in moles
    n = (4/3. * np.pi * r1**3)/(R*T0) * \
        (p0 - 2*mu*(t0/r0)* \
             ((r0/r1) - (r0/r1)**7) * \
             (1 + (1/alfa - 1) * (r1/r0)**2))
    gas_mass = n*M

    # Solve balloon-radius-roots for each height level
    
    # Constants for the polynomials
    a8 = (1/alfa-1)*2*mu*t0/(r0**2)
    a6 = 2*mu*t0
    a2 = -2*mu*t0*(1/alfa-1)*(r0**4)
    a0 = -2*mu*t0*(r0**6)
    
    all_radii = []

    try:
        x, y = p.shape
        
        for i in range(0, x):
            radii = []
            for j in range(0, y):
                a4 = -3*n[j]*R/(4*np.pi)
                # 8th degree polynomial
                poly = [a8,        # r^8
                        p[i,j],    # r^7
                        a6,        # r^6
                        0,         # r^5
                        a4*T[i,j], # r^4
                        0,         # r^3
                        a2,        # r^2
                        0,         # r^1
                        a0]        # r^0

                roots = np.roots(poly)
        
                for r in roots:
                    if r.real > 0 and r.imag == 0:
                        radii.append(r.real)
            all_radii.append(np.array(radii))


    except ValueError:        
        for i in range(0, len(p)):
            a4 = -3*n*R/(4*np.pi)
            # 8th degree polynomial
            poly = [a8,        # r^8
                    p[i],    # r^7
                    a6,        # r^6
                    0,         # r^5
                    a4*T[i], # r^4
                    0,         # r^3
                    a2,        # r^2
                    0,         # r^1
                    a0]        # r^0

            roots = np.roots(poly)
        
            for r in roots:
                if r.real > 0 and r.imag == 0:
                    all_radii.append(r.real)

    all_radii = np.array(all_radii)


    return all_radii, gas_mass


