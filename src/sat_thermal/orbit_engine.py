from math import cos, sin

def _sun_unit_vector(ra_rad, dec_rad):
    """
    Convert spherical sky coordinates to a Cartesian unit vector.
    """
    cos_d = cos(dec_rad)
    return (cos_d * cos(ra_rad),
            cos_d * sin(ra_rad),
            sin(dec_rad))
