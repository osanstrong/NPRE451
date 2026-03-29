import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig

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
            elif "<<DATA>>" in line:
                dat_start = i+1
            elif "<<END>>" in line:
                dat_end = i
        spec_c = np.array([float(l.replace(" ","")) for l in lines[dat_start:dat_end]])
        if not crop is None:
            spec_c = spec_c[:crop]
        return spec_c, time_s
        

def get_spec_cps(path, crop=None):
    spec_c, time_s = get_spec(path, crop=crop)
    return spec_c / time_s


def make_maestro_file(path):
    spec_c, time_s = get_spec(path)
    length = len(spec_c)
    header = f"""$SPEC_ID:
No sample description was entered.
$SPEC_REM:
DET# 1
DETDESC# TL-101B-03 digiBASE-RH SN 16039579
AP# Maestro Version 7.01
$DATE_MEA:
03/06/2026 13:19:43
$MEAS_TIM:
{int(time_s)} 000
$DATA:
"""
    footer = """
$ROI:
$PRESETS:
Real Time
300
0
$ENER_FIT:
0.000000 0.000000
$MCA_CAL:
3
0.000000E+000 0.000000E+000 0.000000E+000 
$SHAPE_CAL:
3
0.000000E+000 0.000000E+000 0.000000E+000"""
    entries = [str(c) for c in spec_c]
    entries[0] = f"{entries[0]}\t{length-1}"
    string = header + "\n".join(entries) + footer
    with open(f"{path}.Spe", "w") as out:
        out.write(string)


def vbar(x, label, c=None, style="dashed"):
    [ymin, ymax] = plt.ylim()
    plt.vlines(x, ymin, ymax, label=label, colors=c, linestyles=style)
    plt.ylim(ymin, ymax)

target_labels = {
    "137Cs":"cesium-137",
    "57Co":"cobalt-57",
    "133Ba":"barium-133",
    "54Mn":"manganese-54",
    "60Co":"cobalt-60",
    "EuUnknown":"europium-152",
    "Background":"background"
}

# Expected E (keV): Centroid, FWHM, Net count rate, error : bf, br
id_peaks = {
    "137Cs":[
        "661.57:1005.71,11.77,93502,476:",
    ],
    "60Co":[
        "1173:1783.94,20.57,21207,412",
        "1332:2026.79,24.59,19193,245"
    ],
    "54Mn":[
        "835:1274.12,4.98,186,28"
    ],
    "133Ba":[
        "81:120.35,3.87,68831,472",
        "276:417.6,6.94,12052,308",
        "302:457.92,7.11,40740,522",
        "356:539.11,7.71,85055,412",
        "384:581.68,8.18,10686,174"
    ],
    "57Co":[
        "136:195.76,1.21,175,27",
        "122:182.01,4.97,1102,41"
    ],
    "EuUnknown":[
        "1408:2143.64,26.3,17031,384:"
    ]
}
half_lives_y = {
    "EuUnknown":13.517
}
calpeaks = [ #Must be 1st in list
    "137Cs",
    "60Co"
    # id_peaks["137Cs"][0],
    # id_peaks["60Co"][0]
]

cy = [float(id_peaks[peak][0].split(":")[0]) for peak in calpeaks]
cx = [float(id_peaks[peak][0].split(":")[1].split(",")[0]) for peak in calpeaks]
m = (cy[1]-cy[0]) / (cx[1]-cx[0])
b = cy[0]-(cx[0]*m)

im = (cx[1] - cx[0]) / (cy[1]-cy[0])
ib = cx[0] - (cy[0]*im)
print(f"Calibration: {m:.5f} x + {b:.5f}")

for target in id_peaks:
    for peak in id_peaks[target]:
        energy = float(peak.split(":")[0])
        bin = im*energy + ib
        print(f"For E = {energy:.5f}keV expect bin {bin:.5f}")


target_specs = {}
target_times = {}
for t in target_labels:
    path = f"dat/{t}.mca"
    spec_c, time_s = get_spec(path)
    target_times[t] = time_s
    target_specs[t] = spec_c / time_s
    # make_maestro_file(path)

# ehhhhhhhhh i'll do the plotting n stuff tomorrow

spec_bins = np.array(range(8192))

test = "137Cs"
# make_maestro_file(f"dat/{test}.mca")

test_spec = target_specs[test]
min_height = 0.02
min_width = 20 #bins
peaks = sig.find_peaks(test_spec, height=min_height, width=min_width)
print(peaks)
plt.plot(spec_bins, target_specs[test], label=target_labels[test])
plt.xlabel("Bin #")
plt.ylabel("Intensity (counts per second)")
plt.legend()
plt.show()