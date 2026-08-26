import numpy as np


def invariant_mass(e1, px1, py1, pz1, e2, px2, py2, pz2):
    e = e1 + e2
    px = px1 + px2
    py = py1 + py2
    pz = pz1 + pz2
    m_squared = e**2 - px**2 - py**2 - pz**2
    return np.sqrt(np.clip(m_squared, 0, None))
