import numpy as np
from flask import Flask, render_template, request, jsonify
from flask_cors import CORS

app = Flask(__name__)
CORS(app, resources={r"/calculate": {"origins": "*"}})


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
    try:
        data = request.get_json()
        h1 = float(data['h1'])
        h2 = float(data['h2'])
        d = float(data['d'])
        freq = float(data['freq'])
        er = float(data['er'])
        pol = data['pol']

        # Call your existing Friis equation function
        result = friis_equation_with_ground_presence(h1, h2, d, freq, er, pol)

        return jsonify(result=result)
    except Exception as e:
        return jsonify(error=str(e)), 500

if __name__ == '__main__':
    app.run(debug=True)
