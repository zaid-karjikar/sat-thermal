from math import cos, sin


def _sun_unit_vector(ra_rad, dec_rad):
    """
    Convert spherical sky coordinates to a Cartesian unit vector.
    """
    cos_d = cos(dec_rad)
    return (cos_d * cos(ra_rad), cos_d * sin(ra_rad), sin(dec_rad))


def _position_eci(
    orbital_radius_m,
    inclination_rad,
    raan_rad,
    true_anomaly_0_rad,
    mean_motion_rad_s,
    t_s,
):
    """
    Returns the ECI position vector (x, y, z) in meters for a circular orbit at time t_s.
    """
    true_anomaly_rad = true_anomaly_0_rad + mean_motion_rad_s * t_s

    # step 1 - position in orbital plane
    x_p = orbital_radius_m * cos(true_anomaly_rad)
    y_p = orbital_radius_m * sin(true_anomaly_rad)
    z_p = 0.0

    # step 2 - tilt by inclination
    cos_i, sin_i = cos(inclination_rad), sin(inclination_rad)
    x_n = x_p
    y_n = y_p * cos_i - z_p * sin_i
    z_n = y_p * sin_i + z_p * cos_i

    # step 3 - rotate by RAAN
    cos_o, sin_o = cos(raan_rad), sin(raan_rad)
    x = x_n * cos_o - y_n * sin_o
    y = x_n * sin_o + y_n * cos_o
    z = z_n

    return (x, y, z)
