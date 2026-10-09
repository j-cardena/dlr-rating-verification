import numpy as np
from dlrverify.ieee738 import ampacity, k_angle


def test_drake_example():
    # IEEE 738-2012 annex: Drake, Ta = 40 C, V = 0.61 m/s perpendicular, Ts = 100 C -> about 1025 A
    G = 22.45 / (0.8 * 0.02814)          # irradiance giving the annex solar gain of ~22.45 W/m
    I = float(ampacity(40.0, 0.61, 90.0, G, Tmax=100.0, elev=0.0))
    assert abs(I - 1025) < 5


def test_angle_factor_limits():
    assert abs(k_angle(90.0) - 1.0) < 0.01      # perpendicular wind
    assert k_angle(0.0) < 0.5                    # parallel wind


def test_monotonic_in_wind_and_temperature():
    I = ampacity(25.0, np.array([0.5, 1.0, 2.0, 5.0]), 90.0, 500.0, Tmax=75.0)
    assert np.all(np.diff(I) > 0)
    J = ampacity(np.array([0.0, 20.0, 35.0]), 2.0, 90.0, 500.0, Tmax=75.0)
    assert np.all(np.diff(J) < 0)
