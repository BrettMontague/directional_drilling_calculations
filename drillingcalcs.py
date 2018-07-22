import math
import numpy as np
from typing import Tuple


def minimum_curvature_calc(
        distance_between_stations: float,
        azimuth_survey_station_1: float,
        dip_survey_station_1: float,
        azimuth_survey_station_2: float,
        dip_survey_station_2: float,
        convergence: float = 0) -> Tuple[float, float, float, float]:
    ''' 
    Using the length between two surveys, an initial azimuth, an initial dip,
    a secondary azimuth, a secondary dip and an optional grid convergence argument
    , return the tvd difference, northing difference, easting difference and 
    dogleg severity in (°/30m).
    See https://www.onepetro.org/journal-paper/SPE-3362-PA for detail on calculations.
    '''

    # Set variables for the inclination and azimuth of the first two survey points
    inclination_1 = math.radians(dip_survey_station_1)
    azimuth_1 = math.radians(azimuth_survey_station_1 - convergence)

    # Set variables for the inclination and azimuth of the second two survey points
    inclination_2 = math.radians(dip_survey_station_2)
    azimuth_2 = math.radians(azimuth_survey_station_2 - convergence)


    # Calculate the DogLeg
    dogleg = np.arccos(
        np.cos(inclination_2 - inclination_1) - np.sin(inclination_1) *
        np.sin(inclination_2) * (1 - np.cos(azimuth_2 - azimuth_1)))

    # Warn if length is zero, prevent divison by zero.
    if distance_between_stations == 0:
        print('''ERROR: Cannot perform calculation on zero length between 
            survey stations. Return values set to 'None' ''')
        return None, None, None, None

    # Calculate the dogleg severity in (°/30m)
    dogleg_severity = math.degrees(dogleg * (30 / distance_between_stations))

    # Calculate ratio factor.  If there's no dogleg, calculate with balanced tangential instead of minimum curvature.
    if dogleg != 0.0:
        ratio_factor = 2 * np.tan(dogleg / 2) / dogleg  # minimum curvature
    else:
        ratio_factor = 1  # balanced tangential

    # Calculation for TVD depth difference between first and second survey stations
    tvd_difference = (0.5 * distance_between_stations * (
        np.cos(inclination_1) + np.cos(inclination_2)) * ratio_factor)

    # Calculation for northing difference between first and second survey stations
    northing_difference = 0.5 * distance_between_stations * (
        np.sin(inclination_1) * np.cos(azimuth_1) +
        np.sin(inclination_2) * np.cos(azimuth_2)) * ratio_factor

    # Calculation for easting difference between first and second survey stations
    easting_difference = 0.5 * distance_between_stations * (
        np.sin(inclination_1) * np.sin(azimuth_1) +
        np.sin(inclination_2) * np.sin(azimuth_2)) * ratio_factor

    return tvd_difference, northing_difference, easting_difference, dogleg_severity


def vertical_section_calc(northing: float, easting: float,
                          vertical_section_difference: float
                          ) -> Tuple[float, float, float, float]:
    ''' 
    using northing, easting and vertical section direction as inputs 
    and return directional difference, closure distance, closure azimuth
    , vertical section difference (vsd). 
    '''

    # Calculation for closure distance
    closure_distance = (northing**2 + easting**2)**.5

    # Calculation for closure azimuth
    if northing == 0:
        closure_azimuth = 0
    else:
        closure_azimuth = math.degrees(np.arctan(easting / northing))

    # Calculation for directional difference
    directional_difference = closure_azimuth - vertical_section_difference

    # Calculation for vertical section
    vertical_section = (closure_distance * (np.cos(
        math.radians(vertical_section_difference - closure_azimuth))))

    return directional_difference, closure_distance, closure_azimuth, vertical_section
