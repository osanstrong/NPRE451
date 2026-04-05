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

DIST = 45 + 5 # mm
DIST_ERR = 1 # mm
DET_RADIUS = 49.1 / 2 # mm
DET_LENGTH = 36 # mm
def average_dist(h, r, c):
    return np.sqrt((r**2)/2 + (h**2)/12 + c**2)
def geom_atten(a, c):
    d = average_dist(DET_LENGTH, DET_RADIUS, c)
    # return a**2/(4*d2)
    return 0.5*(1 - (d/(np.sqrt(d**2 + a**2))))
def geom_error(a, d, d_err):
    Omeg = geom_atten(a, d)
    return Omeg * (2*d_err/d)

GEOM_ATTEN = geom_atten(DET_RADIUS, DIST)
GEOM_ERROR = geom_error(DET_RADIUS, DIST, DIST_ERR)
GEOM_ERROR = 0.00134215660207 # https://www.desmos.com/calculator/t2e4wghni1
print(f"Geometric Efficiency: {GEOM_ATTEN:.5f} +- {GEOM_ERROR:.5f}")



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
        "276:417.56,6.73,10915,144:7.16,100",
        "302:457.92,6.93,26014,198:18.34,100",
        "356:539.11,7.66,82174,316:62.05,100",
        "384:581.65,8.14,10868,121:8.94,100",
    ],
    "57Co":[
        "136:195.76,1.21,175,27:10.68,100",
        "122:182.01,4.97,1102,41:85.60,100"
    ],
    "EuUnknown":[
        # "1408:2143.64,26.3,17031,384:20.87,72.08" # Positron
        # "122:181.84,5.06,181140,1245:", # b+
        # "244:369.22,6.83,32375,968:", # b+
        # "344:521.26,7.72,55067,961:",
        # "779:1183.91,15.09,17175,647:",
        # "964:1466.18,18.52,16109,542:",
        
        # Using data from gammaspectacular
        # "1408:2143.64,26.3,17031,384:21.005,72.08" # b+
        # "122:181.84,5.06,181140,1245:28.58,72.08", # b+
        "244:369.22,6.83,32375,968:7.583,72.08", # b+
        # "344:521.26,7.72,55067,961:26.5,27.92", # b-
        "344:521.26,7.72,76021,492:26.5,27.92", # b-
        # "779:1183.91,15.09,17175,647:12.942,27.92",
        "779:1183.91,15.09,16301,357:12.942,27.92",
        # "964:1466.18,18.52,16109,542:14.605,72.08", 
        # "964:1466.10,18.22,14806,305:14.605,72.08", 
    ]
}
for iso in id_peaks:
    print("______________________________")
    print(f"{iso}:")
    for line in id_peaks[iso]:
        dat_sets = line.split(":")
        N = int(dat_sets[1].split(",")[2])
        Nstr = '{0:,}'.format(N)
        Nerr = '{0:,}'.format(int(dat_sets[1].split(',')[3]))
        E = dat_sets[0]
        print(f"@ {" "*(5-len(E))}{E} kev: N = {" "*(10-len(Nstr))}{Nstr} +- {Nerr}")



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

cy = [float(id_peaks[peak][0].split(":")[0]) for peak in calpeaks]
cx = [float(id_peaks[peak][0].split(":")[1].split(",")[0]) for peak in calpeaks]
m = (cy[1]-cy[0]) / (cx[1]-cx[0])
b = cy[0]-(cx[0]*m)

im = (cx[1] - cx[0]) / (cy[1]-cy[0])
ib = cx[0] - (cy[0]*im)


target_specs = {}
target_times = {}
for t in target_labels:
    path = f"dat/{t}.mca"
    spec_c, time_s = get_spec(path)
    target_times[t] = time_s
    target_specs[t] = spec_c / time_s
    spec_bins = np.array(range(len(spec_c)))
    spec_e = m*spec_bins + b
    # plt.plot(spec_e, spec_c/time_s, label=f"{target_labels[t]}")
    # plt.xlabel("Energy (keV)")
    # plt.ylabel("Intensity (counts per second)")
    # plt.legend()
    # plt.show()
def norm(x, mu, sig):
    return (1/np.sqrt(2*np.pi*sig**2)) * np.exp(-(x-mu)**2 / (2*sig**2))
for iso in id_peaks:
    base_spec = target_specs[iso] - target_specs["Background"]
    for peak in id_peaks[iso]:
        info_sets = peak.split(":")
        E = info_sets[0]

        spec_info = info_sets[1].split(",")
        centroid = float(spec_info[0])
        FWHM = float(spec_info[1])
        N = float(spec_info[2])
        N_err = float(spec_info[3])

        sigma = FWHM / 2.355
        
        range_bins = [int(centroid-3*sigma), int(centroid+7*sigma)+1]
        rel_bins = np.array(range(range_bins[0], range_bins[1]))
        base_peak = base_spec[range_bins[0]:range_bins[1]]
        fit_peak = norm(rel_bins, centroid, sigma)*N / target_times[iso]
        fit_low = fit_peak * (N-N_err)/N
        fit_high = fit_peak * (N+N_err)/N

        rel_E = rel_bins*m + b




print(f"Calibration: {m:.5f} x + {b:.5f}")
for eu_peak in id_peaks["EuUnknown"]:
    centroid = float(eu_peak.split(":")[1].split(",")[0])
    E = m*centroid + b
    print(f"bin {int(centroid)} -> {int(E)} keV")

known_isos = [
    "137Cs", "60Co", 
    "54Mn", 
    "133Ba",
    "57Co",
    ]

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
        # print(E)
        if not E == int(E_int):
            # print(f"{E_int} does not match {E}")
            continue
        # else:
            # print(f"Peak found: {E}")

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
        # print(f"___________________________________________________")
        # print(f"{iso} @ {E} keV:")
        # print(f"A = {BASE_ACTIVITY} * exp(-{lamb_y:.5f} * {elapsed_y:.3f})")
        # print(f"S = {A:.3f} * {T} * ({TB:.3f})")
        # print(f"S = {S:.3f} vs observed N={N:.3f}+-{N_err:.3f}")

        eff = N / S
        err = N_err / S
        return E, eff, err, A
        
    print(f"Unable to find peak info for {iso} {E_int}keV")
    return None, None, None, None

print(f"S = A * T * (BF * BR)")
all_E = []
all_aeff = []
all_aerr = []
all_ieff = []
all_ierr = []

latex_labels={
    "137Cs":"\\CsSrc",
    "60Co":"\\CoSrcA",
    "57Co":"\\CoSrcB",
    "54Mn":"\\MnSrc",
    "133Ba":"\\BaSrc",
}
for iso in known_isos:
    E_list = []
    aeff_list = []
    aerr_list = []
    ieff_list = []
    ierr_list = []
    peaks = id_peaks[iso]
    for peak in peaks:

        info_sets = peak.split(":")
        E = info_sets[0]

        spec_info = info_sets[1].split(",")
        centroid = float(spec_info[0])
        FWHM = float(spec_info[1])
        N = float(spec_info[2])
        N_err = float(spec_info[3])

        branch_info = info_sets[2].split(",")
        BF = float(branch_info[1])/100
        BR = float(branch_info[0])/100 #pretty sure the TBR is tabulated as just branching ratio
        
        jBR = BR / (BF)

        # E = int(peak.split(":")[0])
        E, eff, err, A = abs_eff(iso, E)
        ieff = eff / GEOM_ATTEN
        ierr = ieff * ((err/eff)**2 + (GEOM_ERROR/GEOM_ATTEN)**2)**0.5
        hl_y = half_lives_y[iso]
        # print(f"Absolute efficiency at {E} keV: {eff:.5f} +- {err:.5f}")
        # print(f"Geometric efficiency of {GEOM_ATTEN:.5f} +- {GEOM_ERROR:.5f}")
        # print(f"--> Intrinsic efficiency @ {E} keV: {ieff:.5f} +- {ierr:.5f}")
        A_micCi = A / MICROCURIE
        print(f"{latex_labels[iso]} & {E} & ${int(N)} \\pm {int(N_err)}$ & {BF:.1f} & {jBR:.3f} & {hl_y:.3f} & {A_micCi:.3f} & ${ieff:.4f} \\pm {ierr:.4f}$ \\\\")
        E_list.append(E)
        all_E.append(E)
        aeff_list.append(eff)
        all_aeff.append(eff)
        aerr_list.append(err)
        all_aerr.append(err)
        
        ieff_list.append(ieff)
        all_ieff.append(ieff)
        ierr_list.append(ierr)
        all_ierr.append(ierr)
    E_list = np.array(E_list)
    aeff_list = np.array(aeff_list)*100
    aerr_list = np.array(aerr_list)*100
    ieff_list = np.array(ieff_list)*100
    ierr_list = np.array(ierr_list)*100
    plt.scatter(E_list, ieff_list, label=f"{iso}", marker="x")
# l0, = plt.plot(0, 0, marker='o', color='b')
    (_, caps, _) = plt.errorbar(E_list, ieff_list, yerr=ierr_list, capsize=10, elinewidth=1, fmt="none")
    for cap in caps:
        cap.set_color('black')
        cap.set_markeredgewidth(1)

ln_E = np.log(np.array(all_E))
ln_ieff = np.log(np.array(all_ieff))
all_ierr = np.log(np.array(all_ierr))
N_TERMS = 4
poly, info = np.polynomial.polynomial.Polynomial.fit(ln_E, ln_ieff, N_TERMS-1, domain=[-1,1], w=1/all_ierr, full=True)
print(poly)
print(info)
resid = info[0][0] #/ len(ln_ieff)
ln_ieff_fit_err = resid**0.5
ln_aeff_fit_err = ln_ieff_fit_err * GEOM_ATTEN # Ignore the geometric error because we're discounting that effect entirely
ln_E_range = np.linspace(min(ln_E), max(ln_E), 1000)
ln_ieff_fit = poly(ln_E_range)
print(f"{poly(0)}")
print(f"{poly(2)}")
# print(f"{np.polyval(poly.convert().coef, 0)}")

Eu_Es = []
Eu_effs = []
Eu_As = [] # Activity estimates from each peak
Eu_A_errs = [] # Error for each activity estimate, taking into account the error of the fit and the error of the peak
Eu_ages = []
decays_Eu = []
decay_vars_Eu = []
for peak in id_peaks["EuUnknown"]:
    dat_Eu = peak.split(":")

    E_Eu = int(dat_Eu[0])
    ieff_Eu = np.exp(poly(np.log(E_Eu)))
    ieff_err = ln_ieff_fit_err*ieff_Eu
    aeff_Eu = ieff_Eu * GEOM_ATTEN
    aeff_err = ieff_err * GEOM_ATTEN # Don't factor the geometric error in, we're using the data straight

    Eu_Es.append(E_Eu)
    Eu_effs.append(aeff_Eu)
    # print("_"*30)
    # print(f"Efficiency @ {E_Eu} keV: {aeff_Eu:.5f} +- {aeff_err:.5f}")
    # eff = N / S
    N_Eu = int(dat_Eu[1].split(",")[2])
    Nerr_Eu = int(dat_Eu[1].split(",")[3])
    
    S_Eu = int(N_Eu / aeff_Eu)
    Serr_Eu = S_Eu * ((Nerr_Eu / N_Eu)**2 + (aeff_err/aeff_Eu)**2)**0.5
    # print(f"152Eu @ {E_Eu} keV: N = {N_Eu} +- {Nerr_Eu} -> S = {S_Eu} +- {Serr_Eu}")
    # S = A * T * TBF
    
    
    
    
    BF_Eu = float(dat_Eu[2].split(",")[0]) / 100
    BR_Eu = float(dat_Eu[2].split(",")[1]) / 100
    TB_Eu = BF_Eu
    Time_Eu = 300
    A_Eu = S_Eu / (TB_Eu * Time_Eu)
    A_micCi = A_Eu / MICROCURIE
    A_Eu_err = Serr_Eu / (TB_Eu * Time_Eu)
    A_micCi_err = A_Eu_err / MICROCURIE
    Eu_As.append(A_Eu)
    Eu_A_errs.append(A_Eu_err)
    print(f"{E_Eu} & ${N_Eu} \\pm {Nerr_Eu}$ & ${aeff_Eu:.5f} \\pm {aeff_err:.5f}$ & ${A_micCi:.5f} \\pm {A_micCi_err:.5f}$ \\\\")
    # A_Eu = np.array([S_Eu-Serr_Eu, S_Eu, S_Eu+Serr_Eu]) / (TB_Eu * Time_Eu)
    # decay_Eu = A_Eu / OG_A_Eu
    # decays_Eu.append(decay_Eu[1])
    # decay_vars_Eu.append((decay_Eu[1]-decay_Eu[0])**2)
    # lamb_Eu = np.log(2) / half_lives_y["EuUnknown"]
    # elap_Eu = np.log(decay_Eu) / (-lamb_Eu)
    # Eu_ages.append(elap_Eu)
    # print(f"From A of {OG_A_Eu:,.3f} to {A_Eu} -> {elap_Eu} years elapsed")

net_A = sum(Eu_As) / len(Eu_As)
net_A_err = np.sqrt(sum(err**2 for err in Eu_A_errs)) / len(Eu_A_errs)
print(f"Activity estimate: ${net_A/MICROCURIE:.5f} \\pm {net_A_err/MICROCURIE:.5f}$~\miCi")
OG_A_Eu = 10*MICROCURIE
decay_Eu = net_A / OG_A_Eu
decay_err = net_A_err / OG_A_Eu
lamb_Eu = np.log(2) / half_lives_y["EuUnknown"]
ln_dec = np.log(decay_Eu)
ln_dec_err = decay_err / decay_Eu
elap_Eu = ln_dec / (-lamb_Eu)
elap_err = ln_dec_err / (lamb_Eu)
print(f"Final age estimate: ${elap_Eu:.5f} \\pm {elap_err:.5f}$~years")

Eu_Es = np.array(Eu_Es)
Eu_ln_ieffs_base = poly(np.log(Eu_Es))
Eu_ieffs_lower = np.exp(Eu_ln_ieffs_base - ln_ieff_fit_err)
Eu_ieffs_upper = np.exp(Eu_ln_ieffs_base + ln_ieff_fit_err)

# mean_decay = sum(decays_Eu) / len(decays_Eu)
# var_decay = sum(decay_vars_Eu) / (len(decays_Eu)**2)
# sig_decay = var_decay**0.5
# print(f"sig: {sig_decay}")
# decay_range = np.array([mean_decay+sig_decay, mean_decay, mean_decay-sig_decay])
# lamb_Eu = np.log(2) / half_lives_y["EuUnknown"]
# elap_range = np.log(decay_range) / (-lamb_Eu)
# print(elap_range[1:3]-elap_range[0:2])
# print(f"{elap_range} years passed")

E_range = np.exp(ln_E_range)
ieff_fit = np.exp(ln_ieff_fit)
ieff_lower = np.exp(ln_ieff_fit-ln_ieff_fit_err)
ieff_upper = np.exp(ln_ieff_fit+ln_ieff_fit_err)

plt.plot(np.exp(ln_E_range), np.exp(ln_ieff_fit)*100, label=f"{N_TERMS} term ln fit")
plt.fill_between(E_range, ieff_lower*100, ieff_upper*100, label="Mean Residual", alpha=0.3, zorder=-5)
# plt.plot(np.exp(ln_E), np.exp(poly(ln_E))*100, label="Four term fit 2")
plt.scatter(np.array(Eu_Es), np.array(Eu_effs)*100/GEOM_ATTEN, label="Est. for 152Eu", color=f"C{len(known_isos)}")
# plt.vlines(Eu_Es, Eu_ieffs_lower*100, Eu_ieffs_upper*100, label="Est. for 152Eu", color=f"C{len(known_isos)}")
# (_, caps, _) = plt.errorbar(Eu_Es, np.exp(Eu_ln_ieffs_base), yerr=(Eu_ieffs_upper-Eu_ieffs_lower)/2, capsize=10, elinewidth=1, fmt="none")
# for cap in caps:
#     cap.set_color('black')
#     cap.set_markeredgewidth(1)
plt.ylabel("Intrinsic Efficiency (%)")
plt.xlabel("Energy (keV)")
plt.loglog()
plt.legend()
plt.show()
quit()


for target in id_peaks:
    for peak in id_peaks[target]:
        energy = float(peak.split(":")[0])
        centroid = im*energy + ib
        print(f"For E = {energy:.5f}keV expect bin {centroid:.5f}")




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


DIST = 45 #mm
DIST_S = 1 #also mm, uncertainty
