"""Steady-state ampacity after IEEE Std 738-2012 (SI form), vectorised.

Simplification for the worked example: solar heat gain uses measured global
horizontal irradiance on the projected diameter, q_s = alpha_s * G * D,
instead of the standard's clear-sky model with solar geometry.
"""
import numpy as np

# ACSR 26/7 Drake 795 kcmil, as in the IEEE 738 annex example
DRAKE = dict(D=0.02814, R25=7.283e-5, R75=8.688e-5, eps=0.8, absorp=0.8)


def _air(Tfilm, elev):
    mu = 1.458e-6 * (Tfilm + 273.0) ** 1.5 / (Tfilm + 383.4)
    rho = (1.293 - 1.525e-4 * elev + 6.379e-9 * elev ** 2) / (1 + 0.00367 * Tfilm)
    k = 2.424e-2 + 7.477e-5 * Tfilm - 4.407e-9 * Tfilm ** 2
    return mu, rho, k


def k_angle(phi_deg):
    """Wind direction factor; phi = angle between wind and conductor axis."""
    p = np.radians(phi_deg)
    return 1.194 - np.cos(p) + 0.194 * np.cos(2 * p) + 0.368 * np.sin(2 * p)


def ampacity(Ta, V, phi_deg, G, Tmax=100.0, elev=150.0, c=DRAKE):
    Ta, V, phi_deg, G = map(lambda x: np.asarray(x, float), (Ta, V, phi_deg, G))
    D = c["D"]
    Ts = Tmax
    Tfilm = (Ts + Ta) / 2
    mu, rho, k = _air(Tfilm, elev)
    Re = D * rho * np.maximum(V, 0.0) / mu
    Ka = k_angle(phi_deg)
    dT = Ts - Ta
    qc1 = Ka * (1.01 + 1.35 * Re ** 0.52) * k * dT
    qc2 = Ka * 0.754 * Re ** 0.6 * k * dT
    qcn = 3.645 * rho ** 0.5 * D ** 0.75 * dT ** 1.25
    qc = np.maximum(np.maximum(qc1, qc2), qcn)
    qr = 17.8 * D * c["eps"] * (((Ts + 273) / 100) ** 4 - ((Ta + 273) / 100) ** 4)
    qs = c["absorp"] * np.maximum(G, 0.0) * D
    R = c["R25"] + (c["R75"] - c["R25"]) / 50.0 * (Ts - 25.0)
    return np.sqrt(np.maximum(qc + qr - qs, 0.0) / R)


if __name__ == "__main__":
    # IEEE 738-2012 annex example: Drake, Ta=40 C, V=0.61 m/s perpendicular,
    # Ts=100 C, sea-level-ish; standard reports ~1025 A with its solar model
    # (q_s ~ 22 W/m, i.e. G_eff ~ 980 W/m^2 on projected area).
    I = ampacity(40.0, 0.61, 90.0, 22.45 / (0.8 * 0.02814), Tmax=100.0, elev=0.0)
    print(f"Drake check: {float(I):.0f} A (IEEE 738-2012 example: about 1025 A)")
