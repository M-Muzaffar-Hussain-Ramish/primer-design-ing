from dataclasses import dataclass

@dataclass
class ThermoProfile:
    length: int
    count_A: int
    count_T: int
    count_G: int
    count_C: int
    tm_basic: float
    tm_advanced: float
    gc_percent: float


def thermo_profile(seq: str) -> ThermoProfile:
    s = seq.upper()
    a = s.count("A")
    t = s.count("T")
    g = s.count("G")
    c = s.count("C")
    L = a + t + g + c
    if L == 0:
        return ThermoProfile(0, 0, 0, 0, 0, 0.0, 0.0, 0.0)
    tm_basic = 2 * (a + t) + 4 * (g + c)
    tm_advanced = 64.9 + (41 * ((g + c) - 16.4)) / L
    gc = 100.0 * (g + c) / L
    return ThermoProfile(L, a, t, g, c, float(tm_basic), float(tm_advanced), gc)
