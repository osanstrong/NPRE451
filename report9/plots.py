import matplotlib.pyplot as plt
import numpy as np

# Part 1

# Counting Curve
def vbar(x, label, c=None, style="dashed"):
    [ymin, ymax] = plt.ylim()
    plt.vlines(x, ymin, ymax, label=label, colors=c, linestyles=style)
    plt.ylim(ymin, ymax)
# Plot counts as a function of voltage and select operating voltage
volts = np.array([i*100 for i in range(1, 19)])
counts = np.array([
    0,
    0,
    0,
    1,
    0,
    6,
    7,
    48,
    1894,
    11332,
    22042,
    31192,
    35510,
    35529,
    35473,
    37113,
    61472,
    69100,
])
SELECTED_BIAS = 1400

# plt.plot(volts, counts)
# vbar(SELECTED_BIAS, "Selected Bias: +1400 V", c="C1")
# plt.xlabel("Operating Bias (V)")
# plt.ylabel("Counts Observed over 10 seconds (-)")
# plt.legend()
# plt.show()

# Energy Spectrum
def get_spec(path, crop=None):
     with open(path, 'r', encoding="unicode_escape") as file:
        lines = file.read().splitlines()
        time_s = np.nan
        dat_start = 0
        dat_end = 0
        for i in range(len(lines)):
            line = lines[i]
            if "LIVE_TIME" in line:
                time_s = float(lines[i].split(" ")[-1])
                # print(f"Time: {time_s}s")
            elif "$DATA" in line:
                dat_start = i+2
            elif "$ROI" in line:
                dat_end = i
        spec_c = np.array([float(l.replace(" ","")) for l in lines[dat_start:dat_end]])
        if not crop is None:
            spec_c = spec_c[:crop]
        return spec_c, time_s
        

def get_spec_cps(path, crop=None):
    spec_c, time_s = get_spec(path, crop=crop)
    return spec_c / time_s

spec, time = get_spec("dat/lab9_part1.Spe")
print(f"spec: {spec}")

peaks = np.array([312, 535, 940])

E_He = [0.573, 0.191, 0.764]
E_B = [0.84, 1.47, 0.94*2.310 + 0.06*2.792]
E_Li = [2.73, 2.05, 4.78]

for e in [E_He, E_B, E_Li]:
    print(np.array(e) / peaks)

bins = np.array([i for i in range(len(spec))])
plt.plot(bins, spec, label="PuBe source (thermalized)")
vbar(peaks[0], "E(p), 573 keV", "C1")
vbar(peaks[1], "E(T), 191 keV", "C2")
vbar(peaks[2], "E(sum), 764 keV", "C3")
plt.xlabel("Bin #")
plt.ylabel("Counts")
plt.legend()
plt.show()