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

# Centroid, FWHM, Net count rate
id_peaks = {
    "137Cs":[
        "1005.71,11.77,152.83,"
    ]
}
target_specs = {}
for t in target_labels:
    path = f"dat/{t}.mca"
    target_specs[t] = get_spec_cps(path)
    make_maestro_file(path)

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