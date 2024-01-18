import numpy as np

def friis_equation_with_ground_presence(h1, h2, d, freq, er, pol):
    """
    Calculate the loss of a radio link with ground presence.

    Parameters:
    - h1: Transmitting antenna elevation above ground.
    - h2: Receiving antenna elevation above ground.
    - d: Distance between the two antennas (projected onto the ground plane).
    - freq: Signal frequency in Hz.
    - er: Relative permittivity of the ground.
    - pol: Polarization of the signal ('H' for horizontal, 'V' for vertical).

    Returns:
    - Total_received_energy_dBm: Total received energy in dBm.
    """

    c = 299972458e6
    Gr = 1
    Gt = 1
    Pt = 1e-3
    wavelength = c / freq

    phi = np.arctan((h1 + h2) / d)
    direct_wave = np.sqrt(np.abs(h1 - h2)**2 + d**2)
    refl_wave = np.sqrt(d**2 + (h1 + h2)**2)

    if pol == 'H':
        gamma = (np.sin(phi) - np.sqrt(er - np.cos(phi)**2)) / (np.sin(phi) + np.sqrt(er - np.cos(phi)**2))
    elif pol == 'V':
        gamma = (er * np.sin(phi) - np.sqrt(er - np.cos(phi)**2)) / (er * np.sin(phi) + np.sqrt(er - np.cos(phi)**2))
    else:
        raise ValueError(pol + ' is not a valid polarization')

    length_diff = refl_wave - direct_wave
    cos_phase_diff = np.cos(length_diff * 2 * np.pi / wavelength) * np.sign(gamma)

    Direct_energy = Pt * Gt * Gr * wavelength**2 / ((4 * np.pi * direct_wave)**2)
    reflected_energy = Pt * Gt * Gr * wavelength**2 / ((4 * np.pi * refl_wave)**2) * np.abs(gamma)
    Total_received_energy = Direct_energy + cos_phase_diff * reflected_energy
    Total_received_energy_dBm = 10 * np.log10(Total_received_energy * 1e3)

    return Total_received_energy_dBm

# Example usage:
h1 = 10
h2 = 5
d = 100
freq = 2.4e9
er = 4
pol = 'H'

result = friis_equation_with_ground_presence(h1, h2, d, freq, er, pol)
print(result)
