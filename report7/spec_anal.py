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

# Expected E (keV): Centroid, FWHM, Net count rate, error : bf (energy), br (mode)
id_peaks = {
    "137Cs":[
        "662:1005.71,11.77,93502,476:85.1,100",
    ],
    "60Co":[
        "1173:1783.94,20.57,21207,412:99.85,100",
        "1332:2026.79,24.59,19193,245:99.9826,100"
    ],
    "54Mn":[
        "835:1274.12,4.98,186,28:99.976,100"
    ],
    "133Ba":[
        "81:120.35,3.87,68831,472:32.9,100", # Should we add the 2.65% from the 79.6keV peak i wonder
        "276:417.6,6.94,12052,308:7.16,100",
        "302:457.92,7.11,40740,522:18.34,100",
        "356:539.11,7.71,85055,412:62.05,100",
        "384:581.68,8.18,10686,174:8.94,100",
    ],
    "57Co":[
        "136:195.76,1.21,175,27:10.68,100",
        "122:182.01,4.97,1102,41:85.60,100"
    ],
    "EuUnknown":[
        "1408:2143.64,26.3,17031,384:20.87,72.08" # Positron
    ]
}

BASE_TIME_S = 600 # 1 microcurie samples were taken for 10 minutes
DPY = 365.25
half_lives_y = {
    "137Cs":30.08,
    "60Co":1925.28/DPY,
    "54Mn":312.2/DPY,
    "133Ba":10.551,
    "57Co":271.74/DPY,
    "EuUnknown":13.517,
}

# 3129 for Sep, 3099 for Oct
e_sep = 3129
e_oct = 3099
elapsed_d = {
    "137Cs":e_sep,
    "60Co":e_sep,
    "54Mn":e_oct,
    "133Ba":e_sep,
    "57Co":e_oct,
    "EuUnknown":13.517,
}
calpeaks = [ #Must be 1st in list
    "137Cs",
    "60Co"
    # id_peaks["137Cs"][0],
    # id_peaks["60Co"][0]
]

known_isos = ["137Cs", "60Co", "54Mn", "133Ba", "57Co"]

MICROCURIE = 3.7e4
BASE_ACTIVITY = MICROCURIE
def abs_eff(iso, E_int):
    """
    Returns E, e_abs, uncertainty of e_abs
    """
    peak_infos = id_peaks[iso]
    peak = None
    for p in peak_infos:
        info_sets = p.split(":")
        E = int(info_sets[0])
        if not E == E_int:
            continue

        hl_y = half_lives_y[iso]
        lamb_y = np.log(2) / hl_y
        elapsed_y = elapsed_d[iso] / DPY

        branch_info = info_sets[2].split(",")

        spec_info = info_sets[1].split(",")

        A = BASE_ACTIVITY * np.exp(-lamb_y*elapsed_y)
        T = BASE_TIME_S
        BF = float(branch_info[1]) / 100
        BR = float(branch_info[0]) / 100 #pretty sure the TBR is tabulated as just branching ratio
        TB = BR 

        S = A*T*TB

        N = float(spec_info[2])
        N_err = float(spec_info[3])
        print(f"___________________________________________________")
        print(f"{iso} @ {E} keV:")
        print(f"A = {BASE_ACTIVITY} * exp(-{lamb_y:.5f} * {elapsed_y:.3f})")
        print(f"S = {A:.3f} * {T} * ({TB:.3f})")
        print(f"S = {S:.3f} vs observed N={N:.3f}+-{N_err:.3f}")

        eff = N / S
        err = N_err / S
        return E, eff, err
        
    print(f"Unable to find peak info for {iso} {E_int}keV")
    return None, None, None

print(f"S = A * T * (BF * BR)")
E_list = []
aeff_list = []
aerr_list = []
for iso in known_isos:
    peaks = id_peaks[iso]
    for peak in peaks:
        E = int(peak.split(":")[0])
        E, eff, err = abs_eff(iso, E)
        print(f"Efficiency at {E} keV: {eff:.5f} +- {err:.5f}")
        E_list.append(E)
        aeff_list.append(eff)
        aerr_list.append(err)
E_list = np.array(E_list)
aeff_list = np.array(aeff_list)*100
aerr_list = np.array(aerr_list)*100
plt.scatter(E_list, aeff_list, label="Absolute efficiency", marker="x")
# l0, = plt.plot(0, 0, marker='o', color='b')
(_, caps, _) = plt.errorbar(E_list, aeff_list, yerr=aerr_list, label="Error", capsize=10, elinewidth=1, fmt="none")
for cap in caps:
    cap.set_color('red')
    cap.set_markeredgewidth(1)
plt.ylabel("Absolute Efficiency (%)")
plt.xlabel("Energy (keV)")
plt.legend()
plt.show()
quit()

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
plt.xlabel("Energy (keV)")
plt.ylabel("Intensity (counts per second)")
plt.legend()
plt.show()

DIST = 45 #mm
DIST_S = 1 #also mm, uncertainty
