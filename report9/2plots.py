import matplotlib.pyplot as plt
import numpy as np

def raw_to_spec(path):
    with open(path, "r") as raw:
        return np.array([float(n) for n in raw.read().splitlines()])
    
cs_spec = raw_to_spec("dat/cs137.txt")
cs_bins = np.array([i for i in range(len(cs_spec))])

na_spec = raw_to_spec("dat/na22.txt")
na_bins = np.array([i for i in range(len(na_spec))])


mec2 = 511 # keV, rest mass in units of energy
J_per_MeV = 1.602e-13
J_per_keV = J_per_MeV * 1e-3
keV_per_J = 1/J_per_keV
PI = np.pi

def get_compton_Ee(Eg, theta):
    # Eg in keV
    Ee = Eg - get_compton_Egp(Eg, theta)
    return Ee
def get_compton_Egp(Eg, theta):
    # Eg in keV
    Egp = Eg / (1 + (Eg/mec2)*(1-np.cos(theta)))
    return Egp
def vbar(x, label, c=None, style="dashed"):
    [ymin, ymax] = plt.ylim()
    plt.vlines(x, ymin, ymax, label=label, colors=c, linestyles=style)
    plt.ylim(ymin, ymax)

peaks_E = [
    get_compton_Ee(661.7, PI),
    get_compton_Ee(511, PI)
]
peaks_bins = [
    197,
    134
]
m = (peaks_E[1] - peaks_E[0]) / (peaks_bins[1]-peaks_bins[0])
b = peaks_E[0] - peaks_bins[0]*m

print(f"Calibration: E = {m:.5f}*bin + {b:.5f}")

plt.plot(cs_bins*m+b, cs_spec, label="cesium-137")
vbar(peaks_E[0], label=f"661.7 keV Compton Edge ({peaks_E[0]:.2f} keV)", c="C1")
# plt.plot(na_bins*m+b, na_spec, label="sodium-22")
plt.xlabel("Energy (keV)")
plt.ylabel("Counts")
plt.legend()
plt.show()