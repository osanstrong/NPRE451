import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress as linreg

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

peaks_E = np.array([
    get_compton_Ee(661.7, PI),
    get_compton_Ee(511, PI),
    get_compton_Ee(1200, PI)
])
peaks_bins = np.array([
    197,
    134,
    480
])
# m = (peaks_E[1] - peaks_E[0]) / (peaks_bins[1]-peaks_bins[0])
# b = peaks_E[0] - peaks_bins[0]*m
res = linreg(peaks_bins, peaks_E)
m = res.slope
b = res.intercept

ires = linreg(peaks_E, peaks_bins)
im = ires.slope
ib = ires.intercept

print(f"Calibration: E = {m:.5f}*bin + {b:.5f}, correlation R²={res.rvalue**2}")

# plt.plot(cs_bins*m+b, cs_spec, label="cesium-137")
plt.plot(cs_bins, cs_spec, label="cesium-137")
# vbar(peaks_E[0], label=f"661.7 keV Compton Edge ({peaks_E[0]:.2f} keV)", c="C1")
vbar(peaks_bins[0], label=f"661.7 keV Compton Edge ({peaks_E[0]:.2f} keV)", c="C1")
# plt.xlabel("Energy (keV)")
plt.xlabel("Pulse Integral (channel bin)")
plt.ylabel("Counts")
plt.legend()
# plt.xlim([100, 1150])
plt.show()


# plt.plot(na_bins*m+b, na_spec, label="sodium-22")
plt.plot(na_bins, na_spec, label="sodium-22")
# vbar(peaks_E[1], label=f"511 keV Compton Edge ({peaks_E[1]:.2f} keV)", c="C1")
# vbar(peaks_E[2], label=f"1.2 MeV Compton Edge ({peaks_E[2]:.2f} keV)", c="C2")
vbar(peaks_bins[1], label=f"511 keV Compton Edge ({peaks_E[1]:.2f} keV)", c="C1")
vbar(peaks_bins[2], label=f"1.2 MeV Compton Edge ({peaks_E[2]:.2f} keV)", c="C2")
# plt.xlabel("Energy (keV)")
plt.xlabel("Pulse Integral (channel bin)")
plt.ylabel("Counts")
plt.legend()
# plt.xlim([0, 1500])
plt.show()

plt.scatter(peaks_E[0], peaks_bins[0], label="cesium-137")
plt.scatter(peaks_E[1:], peaks_bins[1:], label="sodium-22")
plt.plot(peaks_E, peaks_E*im + ib, label=f"Linear regression (R²={ires.rvalue**2:.3f})")
plt.xlabel("Energy (keV)")
plt.ylabel("Detector Response (channel)")
plt.legend()
plt.show()

# Experiment 2: spectra through shielding
ex2_paths = [f"./dat/Cf252_80ns/{root}.txt" for root in [
    "bare", "lead", "poly1", "poly3", "poly5"
]]

ex2_names = [
    "No Shielding",
    "Lead (7.71±0.01 mm)",
    "PE (12.8±0.1 mm)",
    "PE (25.6±0.2 mm)",
    "PE (38.4±0.3 mm)"
]

ex2_counts = []
ex2_specs = []
for i in range(len(ex2_paths)):
    path = ex2_paths[i]
    name = ex2_names[i]
    counts = raw_to_spec(path)
    if name == "No Shielding":
        counts*=2
    ex2_specs.append(counts)
    cf = crunch_factor = 1
    crunch_counts = np.array([
        sum(counts[i*cf:i*cf+cf]) for i in range(int(len(counts)/(cf*2)))
    ])

    total_counts = sum(counts)
    print(f"Total counts ({name}): {total_counts:,}")

    ener = np.array(range(len(crunch_counts)))*m*cf + b
    plt.plot(ener, (crunch_counts/max(counts*cf)), label=name)
plt.xlabel("Energy (keV)")
plt.ylabel("Relative counts to maximum")
plt.yscale("log")
plt.legend()
plt.show()

for i in range(len(ex2_paths)):
    name = ex2_names[i]
    spec = ex2_specs[i]
    ener = np.array(range(len(spec)))*m + b
    plt.plot(ener, spec, label=name)
plt.xlabel("Energy (keV)")
plt.ylabel("Counts over 10 minutes")
plt.legend()
plt.show()