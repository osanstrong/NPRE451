import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as stat

SKIP_SPEC = 0 #First 10 entries seem to correspond to negative energy?

def get_spec_cps(path, crop=300):
    with open(path, 'r') as file:
        lines = file.read().splitlines()
        time_s = np.nan
        dat_start = 0
        dat_end = 0
        for i in range(len(lines)):
            line = lines[i]
            if "MEAS_TIM" in line:
                time_raw = lines[i+1].split(" ")
                time_s = float(time_raw[0]) + float(time_raw[1])/(10**len(time_raw[1]))
                # print(f"Time: {time_s}s")
            elif "DATA" in line:
                dat_start = i+2+SKIP_SPEC
            elif "ROI" in line:
                dat_end = i
        spec_c = np.array([float(l.replace(" ","")) for l in lines[dat_start:dat_end]])
        return spec_c[:crop] / time_s
                

def vbar(x, label, c=None, style="dashed"):
    [ymin, ymax] = plt.ylim()
    plt.vlines(x, ymin, ymax, label=label, colors=c, linestyles=style)
    plt.ylim(ymin, ymax)

# Stats
def mean(dat_list) -> float:
    return sum(dat_list) / len(dat_list)

def deviations(dat_list) -> list[float]:
    dat_mean = mean(dat_list)
    return [(n-dat_mean) for n in dat_list]

def sample_variance(dat_list) -> float:
    devs = deviations(dat_list)
    n = len(dat_list)
    return sum([dev*dev for dev in devs]) / (n-1)

def sample_stddev(dat_list) -> float:
    return sample_variance(dat_list)**0.5

def scl(dat_list, factor) -> list:
    return [factor*dat for dat in dat_list]

    
isotopes = ["22Na", "54Mn", "60Co", "133Ba", "137Cs"]
bg_spec = get_spec_cps("dat/Calibration/background(15MIN).Spe")
bg_chan = np.array([i for i in range(len(bg_spec))])
iso_specs = {}
for iso in isotopes:
    iso_specs[iso] = get_spec_cps(f"dat/Calibration/{iso}calib.Spe") - bg_spec


iso_cals = { # keV / channel
    # "22Na-g": 511 / iso_specs["22Na"].argmax(),
    "22Na" : 1274.5 / (173-SKIP_SPEC),
    "54Mn" : 834.8 / (115-SKIP_SPEC),
    "60Co" : 1173.2 / (iso_specs["60Co"].argmax()+1),
    "133Ba" : 356 / (50.5-SKIP_SPEC),
    "137Cs" : 661.657 / (iso_specs["137Cs"].argmax()+1)
}

calraws = [iso_cals[key] for key in iso_cals]
calmean = mean(calraws)
calsdev = sample_stddev(calraws)

# print(iso_cals)
print(f"Calibration coeff: {calmean}+-{calsdev} keV/bin")
# print(bg_spec)
# plt.plot(bg_chan, bg_spec, label="Background")
# for iso in isotopes:
#     iso_spec = get_spec_cps(f"dat/Calibration/{iso}calib.Spe")
#     plt.plot(bg_chan, iso_spec, label=f"{iso}")
# plt.xlabel("Channel Number")
# plt.ylabel("Counts per second")
# plt.legend()
# plt.show()

bg_chan = np.array([i for i in range(len(bg_spec))])
# # plt.plot(bg_chan, bg_spec, label="Background")
# for iso in isotopes:
#     iso_spec = iso_specs[iso]
#     plt.plot(bg_chan, iso_spec, label=f"{iso}")

# plt.plot(bg_chan, iso_specs["54Mn"], label="Mn-54")
# plt.plot(bg_chan, iso_specs["137Cs"], label="Cs-137")
# plt.xlabel("Channel Number")
# plt.ylabel("Counts per second")
# plt.legend()
# plt.show()


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

def singE(Egp):
    return Egp - mec2
def doubE(Egp):
    return Egp - 2*mec2

def get_FWHM_I_cent(spec:np.ndarray, peak_idx, cal, radius=100):
    ifrom = max(0, int(peak_idx) - radius)
    ito = int(peak_idx) + radius
    print(f"from: {ifrom}, to: {ito}, spec of len {len(spec)}")
    cent = spec[ifrom:ito].argmax() + ifrom
    I = spec[ifrom:ito].max()
    hm = I*0.5
    left = cent
    right = cent
    for i in range(ifrom):
        left = cent - i
        if spec[left] < hm:
            break

    for i in range(ito, len(spec)):
        right = i
        if spec[right] < hm:
            break

    if left==cent or right==cent:
        print(f"Couldn't find halfsides for E={peak_idx*cal} keV within radius {radius*cal} keV")
    return (right - left)*cal, sum(spec[left:right])*cal, cent*cal

def fwhm_to_sigma(FWHM):
    # FWHM = 2sqrt(2ln2) sig
    return FWHM / (2*np.sqrt(2*np.log(2)))



mph = 60
mpd = mph*24
mpy = mpd*365.25

dpy = 365.25
iso_hlm = {
    "137Cs":30.17*mpy,
    "22Na":2.6*mpy,
    "54Mn":312*mpd,
    "60Co":5.27*mpy,
    "133Ba":10.5*mpy,
}
iso_hly = {
    "137Cs":30.17,
    "22Na":2.6,
    "54Mn":312/dpy,
    "60Co":5.27,
    "133Ba":10.5,
}
# Full Energy Gamma
# "Backscatter Peak": get_compton_Egp(, PI)
# "Compton Edge": get_compton_Ee(, PI), 
# "Single-escape peak": mec2,
# "Double-escape peak": 2*mec2
iso_specs["Unknown"] = get_spec_cps("dat/Unknown/Unknownspectrum.Spe")
iso_cals["Unknown"] = calmean
plot_tomark = {
    # "137Cs":"fcb",
    # "22Na":["fbcam", "c"],
    # "60Co":["fcb", "fbcs"],
    # "54Mn":"fc",
    # "133Ba":["fb"] + ["f"]*3,
    # "Unknown":"fcb"
}
iso_E = {
    "137Cs":661.657,
    "22Na":[1274.5, 511],
    "60Co":[1173.2, 1332.5],
    "54Mn":834.8,
    "133Ba":[356.0, 276, 302, 382],
    "Unknown":661
}
iso_tex = {
    "137Cs":"\\CsSrc",
    "22Na":"\\NaSrc",
    "60Co":"\\CoSrc",
    "54Mn":"\\MnSrc",
    "133Ba":"\\BaSrc",
    "Unknown":"Unknown"
}
# Taken directly from graphs, so need to subtract skip
iso_fwhm = { # bins
    "137Cs":7.93,
    "22Na":[11.26, np.nan],
    "60Co":[9.33, 9.77],
    "54Mn":7.94,
    "133Ba":[5.85, 2.61, 4.00, 3.16, np.nan],
    "Unknown":7.70
}
iso_I = { # cps
    "137Cs":1901.15,
    "22Na":[105.08, np.nan],
    "60Co":[194.66, 163.59],
    "54Mn":2.37,
    "133Ba":[1066.83, 1319.28, 58.74, 60.61, np.nan],
    "Unknown":850.89
}
iso_cent = { # bins
    "137Cs":92.11,
    "22Na":[172.89, np.nan],
    "60Co":[158.19, 179.67],
    "54Mn":114.76,
    "133Ba":[50.96, 41.43, 41.43, np.nan],
    "Unknown":91.12
}

cal_E = {
    "137Cs":661.657,
    "22Na":[1274.5, 511],
    "60Co":[1173.2, 1332.5],
    "54Mn":834.8,
    "133Ba":[356.0, 30.85, 81, 276, 302],
    # "Unknown":661
}
cal_cent = { # bins
    "137Cs":92.11,
    "22Na":[172.89, 71.60],
    "60Co":[158.19, 179.67],
    "54Mn":114.76,
    "133Ba":[50.96, 13.51, 17.90, 41.43, 41.43],
    # "Unknown":91.12
}

# Calibration round 2: linear regression]
cal_El = []
cal_cl = []
for iso in cal_E:
    E = cal_E[iso]
    cent = cal_cent[iso]
    
    if type(E) is not list:
        E = [E]
        cent = [cent]
    # for Ei in E: 
    for i in range(len(E)):
        Ei = E[i]
        ci = (cent[i] - SKIP_SPEC)
        cal_El.append(Ei)
        cal_cl.append(ci)

cal_El = np.array(cal_El)
cal_cl = np.array(cal_cl)
cal_reg = stat.linregress(cal_cl, cal_El)
mcal = cal_reg.slope
bcal = cal_reg.intercept
r2cal = cal_reg.rvalue**2

# plt.scatter(cal_cl, cal_El, label="Identifiable Peaks")
# for iso in cal_E:
#     plt.scatter(cal_cent[iso], cal_E[iso], label=f"{iso}, identified")
# plt.plot(cal_cl, cal_cl*mcal + bcal, label=f"Linear Fit (m={mcal:.3f}, b={bcal:.3f}, r²={r2cal:.3f})")
# plt.xlabel("Peak Centroid (bin)")
# plt.ylabel("Expected Peak Energy (keV)")
# plt.legend()
# plt.show()

spec_Ecal = bg_chan*mcal + bcal
print(f"Calibrated: m={mcal:.3f}, b={bcal:.3f}")
plot_tomark = {
    # "137Cs":"fcb",
    # "22Na":["fbcam", "c"],
    # "60Co":["fcb", "fbcs"],
    # "54Mn":"fc",
    # "133Ba":["fb"] + ["f"]*4,
    # "Unknown":"fcb"
    "Background":["f"]*2
}
iso_specs["Background"] = bg_spec
cal_E["Background"] = np.array([90.43, 197.96]) * mcal + bcal
cal_E["Unknown"] = 661.657

plot_markers = {}
for iso in plot_tomark:
    m = plot_tomark[iso]
    E = cal_E[iso]
    marks = {}
    plot_markers[iso] = marks
    if type(E) is list or type(E) is np.ndarray:
        El = E
        ml = m
        for i in range(len(El)):
            E = El[i]
            m = ml[i]
            plot_markers[iso] = marks
            if iso == "Background":
                if 'f' in m: marks[f"Full-energy gamma ({E:.2f}+-7.69 keV)"] = E
            else:
                if 'f' in m: marks[f"{E} keV Full-energy gamma"] = E
            if 'b' in m: marks[f"{E} keV Backscatter"] = get_compton_Egp(E, PI)
            if 'c' in m: marks[f"{E} keV Compton edge"] = get_compton_Ee(E, PI)
            if 's' in m: marks[f"{E} keV Single escape"] = E - mec2
            if 'd' in m: marks[f"{E} keV Double escape"] = E - 2*mec2
            if 'a' in m: marks["Annihilation (511 keV)"] = mec2
            if 'm' in m: marks[f"{E} keV Sum ({E+mec2} keV)"] = E+mec2
    else:
        if iso == "Unknown":
            if 'f' in m: marks[f"Full-energy gamma (654.67+-7.69 keV)"] = E
        elif iso == "Background":
            if 'f' in m: marks[f"Full-energy gamma ({E:.2f}+-7.69 keV)"] = E
        else:
            if 'f' in m: marks[f"Full-energy gamma ({E} keV)"] = E
        if 'b' in m: marks["Backscatter"] = get_compton_Egp(E, PI)
        if 'c' in m: marks["Compton edge"] = get_compton_Ee(E, PI)
        if 's' in m: marks["Single escape"] = E - mec2
        if 'd' in m: marks["Double escape"] = E - 2*mec2
        if 'a' in m: marks["Annihilation (511 keV)"] = mec2

for target in plot_tomark:
    spec = iso_specs[target]
    pois = plot_markers[target]
    ci = 0

    plt.plot(spec_Ecal, spec, label=target, c=f"C{ci}")
    # plt.fill_betweenx(spec, bg_chan*(calmean-calsdev), bg_chan*(calmean+calsdev), label="Calibration Uncertainty", alpha=0.4)
    for poi in pois:
        ci+=1
        vbar(pois[poi], label=poi, c=f"C{ci}")
    plt.xlabel("Channel Number")
    plt.ylabel("Counts per second")
    plt.legend()
    plt.show()


Rlist = []
Elist = []
El2 = []
clist = []
Ilist = []
time_y = 9 # 9 years from 2017 to 2026
for iso in iso_E:
    # cal = iso_cals[iso]
    err = mcal
    E = iso_E[iso]
    fwhm = iso_fwhm[iso]
    I = iso_fwhm[iso]
    cent = iso_cent[iso]
    
    if type(E) is not list:
        E = [E]
        fwhm = [iso_fwhm[iso]]
        I = [iso_I[iso]]
        cent = [iso_cent[iso]]
    # for Ei in E: 
    for i in range(len(E)):
        Ei = E[i]
        fi = fwhm[i] * mcal
        ci = (cent[i] - SKIP_SPEC) * mcal + bcal
        Ii = I[i] # Don't multiply by bin width bc it's a histogram not direct function of energy
        if Ei == 511:
            continue # Not a photopeak
        if not iso=="Unknown":
            lamb_y = np.log(2) / iso_hly[iso]
            dec_corr = np.exp( time_y * lamb_y)
            # print(f"Corrected by {dec_corr}")
            Ii *= dec_corr
            El2.append(Ei)
            clist.append(ci)
            Ilist.append(Ii)
            # print(f"{iso}, {Ei}: {Ii}")
        # else: 
        #     continue
        # if not Ei == 511:
            # fwhm, I, cent = get_FWHM_I_cent(iso_specs[iso], Ei/cal, cal)
        sigma = fwhm_to_sigma(fi)
        if not np.isnan(fi):
            print(f"{iso_tex[iso]} & {Ei}~keV & ${ci:.2f}\\pm{mcal:.2f}$~keV & {fi:.2f}~keV & {sigma:.2f} \\\\")
            R = fi / ci
            Rlist.append(R)
            Elist.append(ci)
        else:
            print(f"{iso_tex[iso]} & {Ei}~keV & \\na & \\na & \\na \\\\")
# plt.scatter(El2, Ilist, label="Observed Photopeaks")
# for iso in iso_E:
#     E = iso_E[iso]
#     I = iso_fwhm[iso]
#     cent = iso_cent[iso]
    
#     if type(E) is not list:
#         E = [E]
#         I = [iso_I[iso]]
#         cent = [iso_cent[iso]]
#     I_corr = []
#     cent_corr = []
#     # for Ei in E: 
#     for i in range(len(E)):
#         Ei = E[i]
#         ci = (cent[i] - SKIP_SPEC) * mcal + bcal
#         Ii = I[i] # Don't multiply by bin width bc it's a histogram not direct function of energy
#         if Ei == 511:
#             continue # Not a photopeak
#         if not iso=="Unknown":
#             lamb_y = np.log(2) / iso_hly[iso]
#             dec_corr = np.exp( time_y * lamb_y)
#             # print(f"Corrected by {dec_corr}")
#             Ii *= dec_corr
#             I_corr.append(Ii)
#             cent_corr.append(ci)
#     if not len(I_corr) == 0: plt.scatter(cent_corr, I_corr, label=f"{iso}")
# plt.xlabel("Expected Energy (keV)")
# plt.ylabel("Intensity (Counts per second, decay corrected to 2017)")
# plt.yscale("log")
# plt.legend()
# plt.show()

# Elist = np.array(Elist)
# Rlist = np.array(Rlist)
# print(Elist)
# print(Rlist)
# res = stat.linregress(np.log(Elist), np.log(Rlist))
# m = res.slope
# b = res.intercept
# x = np.array([min(np.log(Elist)), max(np.log(Elist))])
# plt.scatter(np.log(np.array(Elist)), np.log(Rlist), label="Observed Photopeaks")
# plt.plot(x, m*x+b, label=f"Linear Regression, Slope: {m:3f}, r²={res.rvalue**2:3f}", c="C1")
# plt.xlabel("Log Observed Centroid Energy (ln(keV))")
# plt.ylabel("Log Observed Resolution (ln(keV))")
# plt.legend()
# plt.show()




# plot_markers = {}
# for iso in plot_tomark:
#     m = plot_tomark[iso]
#     E = iso_E[iso]
#     marks = {}
#     plot_markers[iso] = marks
#     if type(E) is list:
#         El = E
#         ml = m
#         for i in range(len(El)):
#             E = El[i]
#             m = ml[i]
#             plot_markers[iso] = marks
#             if 'f' in m: marks[f"{E} keV Full-energy gamma"] = E
#             if 'b' in m: marks[f"{E} keV Backscatter"] = get_compton_Egp(E, PI)
#             if 'c' in m: marks[f"{E} keV Compton edge"] = get_compton_Ee(E, PI)
#             if 's' in m: marks[f"{E} keV Single escape"] = E - mec2
#             if 'd' in m: marks[f"{E} keV Double escape"] = E - 2*mec2
#             if 'a' in m: marks["Annihilation (511 keV)"] = mec2
#             if 'm' in m: marks[f"{E} keV Sum ({E+mec2} keV)"] = E+mec2
#     else:
#         if iso == "Unknown":
#             if 'f' in m: marks[f"Full-energy gamma (645-663 keV)"] = E
#         else:
#             if 'f' in m: marks[f"Full-energy gamma ({E} keV)"] = E
#         if 'b' in m: marks["Backscatter"] = get_compton_Egp(E, PI)
#         if 'c' in m: marks["Compton edge"] = get_compton_Ee(E, PI)
#         if 's' in m: marks["Single escape"] = E - mec2
#         if 'd' in m: marks["Double escape"] = E - 2*mec2
#         if 'a' in m: marks["Annihilation (511 keV)"] = mec2

# plot_markers_old = {
#     "137Cs":{
#         "Backscatter Peak": get_compton_Egp(661.657, PI),
#         "Compton Edge": get_compton_Ee(661.657, PI),
#         "Full-energy gamma Peak": 661.657,
#         # "Single-escape peak": mec2,
#         # "Double-escape peak": 2*mec2,
#     },
#     "22Na":{
#         "Full-energy gamma Peak": 1274.5,
#         "Backscatter Peak": get_compton_Egp(1274.5, PI),
#         "Compton Edge": get_compton_Ee(1274.5, PI), 
#         "Single-escape Peak": singE(1274.5),
#         "Double-escape Peak": doubE(1274.5),
#         "Annihilation Peak": mec2,
#     },
#     "60Co":{
#         "Full-energy gamma Peak": 1173.2,
#         "Backscatter Peak": get_compton_Egp(1173.2, PI),
#         "Compton Edge": get_compton_Ee(1173.2, PI), 
#         "Single-escape Peak": mec2,
#         "Double-escape Peak": 2*mec2
#     },
#     "54Mn":{
#         "Full-energy gamma Peak": 834.8,
#         "Backscatter Peak": get_compton_Egp(834.8, PI),
#         "Compton Edge": get_compton_Ee(834.8, PI), 
#         # "Single-escape Peak": mec2,
#         # "Double-escape Peak": 2*mec2
#     },
#     "133Ba":{
#         "Full-energy gamma Peak": 356,
#         "Backscatter Peak": get_compton_Egp(356, PI),
#         "Compton Edge": get_compton_Ee(356, PI), 
#         # "Single-escape Peak": mec2,
#         # "Double-escape Peak": 2*mec2
#     },
# }

# # target = "60Co"
# for target in plot_tomark:
#     spec = iso_specs[target]
#     pois = plot_markers[target]
#     ci = 0

#     plt.plot(bg_chan*iso_cals[target], spec, label=target, c=f"C{ci}")
#     # plt.fill_betweenx(spec, bg_chan*(calmean-calsdev), bg_chan*(calmean+calsdev), label="Calibration Uncertainty", alpha=0.4)
#     for poi in pois:
#         ci+=1
#         vbar(pois[poi], label=poi, c=f"C{ci}")
#     plt.xlabel("Channel Number")
#     plt.ylabel("Counts per second")
#     plt.legend()
#     plt.show()
