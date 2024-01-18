import numpy as np

from flask import Flask, render_template, request

app = Flask(__name__)

def friis_equation_with_ground_presence(h1, h2, d, freq, er, pol):
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

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/calculate', methods=['POST'])
def calculate():
    h1 = float(request.form['h1'])
    h2 = float(request.form['h2'])
    d = float(request.form['d'])
    freq = float(request.form['freq'])
    er = float(request.form['er'])
    pol = request.form['pol']

    result = friis_equation_with_ground_presence(h1, h2, d, freq, er, pol)
    return str(result)

if __name__ == '__main__':
    app.run(debug=True)
