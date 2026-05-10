import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

chpath = lambda channel_num: f"dat_72hr/ch{channel_num}.CSV"
ch_df = lambda channel_num: pd.read_csv(chpath(channel_num), sep=";")
timeset = lambda channel_num: np.array(ch_df(channel_num)[0])
def get_realtimes(channel_num):
    


    path = f"dat_72hr/ch{channel_num}_TIMETAG.npy"
    if not os.path.exists(path):
        # df = ch_df(channel_num)

        raw_path = chpath(channel_num)
        # cols = pd.read_csv(raw_path, sep=';', nrows=1)
        dat = pd.read_csv(raw_path, sep=';', skiprows=1)
        ts = np.array(dat.iloc[:,2])
        # print(df.columns)
        # ts = timeset(channel_num)

        np.save(path, ts)
    else:
        # print("getting from file")
        ts = np.load(path)
    return ts[~np.isnan(ts)]
ch0 = get_realtimes(0) / 1000 # ns
ch3 = get_realtimes(3) / 1000
ch4 = get_realtimes(4) / 1000
np.set_printoptions(linewidth=400, edgeitems=5)
print(ch0)
print(ch3)
print(ch4)

ch0 = list(ch0)
ch3 = list(ch3)
ch4 = list(ch4)

ch3_dt = []
ch4_dt = []
j = 0
k = 0

GATE_NS = 600
CH4_0_OFFSET = -11.25 # ns

for i in range(len(ch0)):
    t0 = ch0[i]

    if j < len(ch3):
        t3 = ch3[j]
        dt3 = t3-t0
        while dt3 < -GATE_NS: # This would mean incorrectly interpreting something
            j+=1
            t3 = ch3[j]
            dt3 = t3-t0
            print(f"invalid dt of {dt3} found in channel 3")
        if -GATE_NS < dt3 and dt3 < GATE_NS:
            ch3_dt.append(dt3)
            j+=1

    if k < len(ch4):
        t4 = ch4[k]
        dt4 = t4-t0
        while dt4 < -GATE_NS: # This would mean incorrectly interpreting something
            k+=1
            t4 = ch4[k]
            dt4 = t4-t0
            print(f"invalid dt of {dt4} found in channel 4")
        if -GATE_NS < dt4 and dt4 < GATE_NS:
            ch4_dt.append(dt4)
            k+=1


ch3_dt = np.array(ch3_dt) # ns
ch4_dt = np.array(ch4_dt) + CH4_0_OFFSET

# E = 1/2 m v2
NEUTRON_MASS = 1.675e-27 # kg
EV_PER_J = 6.242e+18
L = 1 # m
S_PER_NS = 1e-9

def dt_to_E(dt_spec_ns):
    v = L / (dt_spec_ns * S_PER_NS)
    E = 0.5 * NEUTRON_MASS * (v**2) * EV_PER_J
    return E / 1000 # keV

ch3_E = dt_to_E(ch3_dt)
ch4_E = dt_to_E(ch4_dt)
shielded_E = ch4_E

# Load existing dt spec
def load_barespec(channel_num):
    path = f"dat_Bare/ch{channel_num}_TOFSPEC.txt"
    with open(path, "r") as tofspec:
        text = tofspec.read()
        vals = [int(s) for s in text.split("\n") if not len(s)==0]
        # print(vals)
        TOF_MAX = 1024
        bins = np.linspace(-TOF_MAX, TOF_MAX, int(len(vals)))
        # bins = np.linspace(-1, 1, len(vals)) # start with just a scale of 0 to 1 then rescale after building a new series of mock data
        # mock
        # pos_idx = np.where(bins > 0) and np.where(bins < GATE_NS)

        return bins, np.array(vals)#[int(len(vals)/2):]


scale_tofspec = lambda dt: dt*1042+44
def load_bare_mockdt(channel_num):
    '''T spectra are loaded in histogram/spectrum form already,
    so to rebin in coarser bins (this doesn't work for finer bins),
    we simply add one 'value' of halfway between the two bins to a new list of mock data'''
    path = f"dat_Bare/ch{channel_num}_TOFSPEC.txt"
    with open(path, "r") as tofspec:
        text = tofspec.read()
        counts = [int(s) for s in text.split("\n") if not len(s)==0]
        bins = np.linspace(-1, 1, len(counts)+1) # start with just a scale of 0 to 1 then rescale after building a new series of mock data
        db = bins[1] - bins[0]
        vals = [] # mock data
        for i in range(len(counts)):
            vals.extend([bins[i]+0.5*db] * counts[i])
        
        return scale_tofspec(np.array(vals))


def load_22hrShield_mockdt(channel_num):
    '''T spectra are loaded in histogram/spectrum form already,
    so to rebin in coarser bins (this doesn't work for finer bins),
    we simply add one 'value' of halfway between the two bins to a new list of mock data'''
    path = f"dat_22hr/ch{channel_num}_TOFSPEC.txt"
    with open(path, "r") as tofspec:
        text = tofspec.read()
        counts = [int(s) for s in text.split("\n") if not len(s)==0]
        bins = np.linspace(-1, 1, len(counts)+1) # start with just a scale of 0 to 1 then rescale after building a new series of mock data
        db = bins[1] - bins[0]
        vals = [] # mock data
        for i in range(len(counts)):
            vals.extend([bins[i]+0.5*db] * counts[i])
        
        return scale_tofspec(np.array(vals))


# bare_dtbins, bare_counts = load_barespec(3)
# bare_Ebins = np.flip(dt_to_E(bare_dtbins))
# print(bare_Ebins)
# bare_counts = np.flip(bare_counts)
# print(f"Max count: {max(bare_counts)}")


bare_dt4 = load_bare_mockdt(4) + CH4_0_OFFSET
bare_dt3 = load_bare_mockdt(3)
bare_E3 = dt_to_E(bare_dt3)
bare_E4 = dt_to_E(bare_dt4)

sh22_dt4 = load_22hrShield_mockdt(4) + CH4_0_OFFSET
sh22_dt3 = load_22hrShield_mockdt(3)
sh22_E4 = dt_to_E(sh22_dt4)
sh22_E3 = dt_to_E(sh22_dt3)



print(f"Lengths: 0: {len(ch0)}, 3: {len(ch3)}, 4: {len(ch4)}")

MIN_E = 15 # 15 keV
MAX_E = 35 # 35 keV

LEFT_E = 20 
RIGHT_E = 30

MIN_E = 0
MAX_E = 200

DT_PAD = 30

LEFT_DT = 417.2 # ns
RIGHT_DT = 511.2

MIN_DT = LEFT_DT - DT_PAD*2
MAX_DT = RIGHT_DT + DT_PAD*0.5

MIN_DT = -GATE_NS
MAX_DT = GATE_NS


UC_DIST = 0.01 # m
UC_TIME = 7 # ns
res = lambda dt: 2*np.sqrt((UC_DIST/L)**2 + (UC_TIME/dt)**2)

res_30 = res(LEFT_DT)
res_20 = res(RIGHT_DT)

UC_E30 = res(LEFT_DT) * 30 / 2
UC_E20 = res(RIGHT_DT) * 20 / 2

print(f"Energy Resolution: {res_20*100:.2f}% to {res_30*100:.2f}%")
print(f"Energy UC: {UC_E20:.2f} keV to {UC_E30:.2f} keV")


NUM_BINS = 50
def hist(dat, label, minx, maxx, scale:float=1, style="solid"):
    bins = np.linspace(minx, maxx, NUM_BINS)
    counts, bins = np.histogram(dat, bins)
    plt.stairs(counts*scale, bins, label=label, linestyle=style)
    return counts, bins

def histE(dat, label, scale:float=1, style="solid"):
    return hist(dat, label, MIN_E, MAX_E, scale=scale, style=style)

def histDT(dat, label, scale:float=1, style="solid"):
    return hist(dat, label, MIN_DT, MAX_DT, scale=scale, style=style)

def vbar(x, label, c=None, style="dashed", err=0, alpha=0.3):
    ymin, ymax = plt.ylim()
    if err == 0:
        plt.plot([x, x], [ymin, ymax], label=label, c=c, linestyle=style)
    else:
        plt.fill_betweenx([ymin, ymax], [x-err]*2, [x+err]*2, label=label, color=c, linestyle=style, alpha=alpha)
    plt.ylim(ymin, ymax)
# hist(ch4, label="Channel 4")

BASE_TIME = 72 # hours
actual_time_fraction = 301 / 2000 # Based on raw data sizes of Ch3, which was 2000MB for the known 72 hours vs 301 for the one that got cut short
SHORT_TIME = BASE_TIME * actual_time_fraction
SHORT_TIME = BASE_TIME * len(ch3_dt) / len(bare_dt3)

# histDT(ch4_dt, label=f"37.5 mm ({SHORT_TIME:.2f} hr), Ch 4", scale=1/SHORT_TIME)
# histDT(ch3_dt, label=f"37.5 mm ({SHORT_TIME:.2f} hr), Ch 3", scale=1/SHORT_TIME)


CUT_TIME = 22.5 # hours
CUT_TIME = BASE_TIME * len(sh22_dt3) / len(bare_dt3)
print(f"Cut time estimate: {CUT_TIME:.3f} h")
# histDT(sh22_dt3, label=f"37.5mm ({CUT_TIME:.2f} hr), Ch3", scale=1/CUT_TIME)
# histDT(sh22_dt4, label=f"37.5mm ({CUT_TIME:.2f} hr), Ch4", scale=1/CUT_TIME)



sh28_dt3 = np.concatenate([ch3_dt, sh22_dt3])
sh28_dt4 = np.concatenate([ch4_dt, sh22_dt4])
sh28_E3 = dt_to_E(sh28_dt3)
sh28_E4 = dt_to_E(sh28_dt4)

#### dt bins
# histDT(sh28_dt3, label=f"37.5mm ({(CUT_TIME+SHORT_TIME):.2f} hr), Ch3", scale=1/(CUT_TIME+SHORT_TIME))
# sh284_counts, bins = histDT(sh28_dt4, label=f"37.5mm ({(CUT_TIME+SHORT_TIME):.2f} hr), Ch4", scale=1/(CUT_TIME+SHORT_TIME))
# histDT(bare_dt3, label="Bare (72 hr), Ch 3", scale=1/BASE_TIME)
# bare4_counts, bins = histDT(bare_dt4, label="Bare (72 hr), Ch 4", scale=1/BASE_TIME)
# # vbar(RIGHT_DT, f"{RIGHT_DT} ns ({RIGHT_E} keV)")
# # vbar(LEFT_DT, f"{LEFT_DT} ns ({LEFT_E} keV)")

# plt.xlabel("dt (ns)")
# plt.ylabel("Counts / hour")
# plt.legend()
# plt.show()

# plt.stairs(sh284_counts / bare4_counts, bins, label=" 37.5mm Fe / Bare (CLYC, Ch 4)")
# vbar(RIGHT_DT, f"{RIGHT_DT} ns ({RIGHT_E} keV)")
# vbar(LEFT_DT, f"{LEFT_DT} ns ({LEFT_E} keV)")
# plt.xlabel("dt (ns)")
# plt.ylabel("Ratio of count rates")
# plt.legend()
# plt.show()

#### E bins


# histE(ch3_E, label=f"Bare Stilbene ({(SHORT_TIME):.2f} hr), Ch3", scale=1/(SHORT_TIME))
# histE(ch4_E, label=f"37.5mm Fe CLYC ({(SHORT_TIME):.2f} hr), Ch4", scale=1/(SHORT_TIME))
# histE(sh22_E3, label=f"Bare Stilbene ({(CUT_TIME):.2f} hr), Ch3", scale=1/(CUT_TIME))
# histE(sh22_E4, label=f"37.5mm Fe CLYC ({(CUT_TIME):.2f} hr), Ch4", scale=1/(CUT_TIME))

MIN_E = 15
MAX_E = 35

histE(sh28_E3, label=f"Bare Stilbene (Total {(CUT_TIME+SHORT_TIME):.2f} hr), Ch3", scale=1/(CUT_TIME+SHORT_TIME))
sh284_counts, bins = histE(sh28_E4, label=f"37.5mm Fe CLYC (Total {(CUT_TIME+SHORT_TIME):.2f} hr), Ch4", scale=1/(CUT_TIME+SHORT_TIME))

histE(bare_E3, label="Bare Stilbene (72 hr), Ch 3", scale=1/BASE_TIME)
bare4_counts, bins = histE(bare_E4, label="Bare CLYC (72 hr), Ch 4", scale=1/BASE_TIME)
vbar(LEFT_E, f"{LEFT_E} keV", err=UC_E20, c="C4")
vbar(25, "25 keV", err=res(0.5*(LEFT_DT+RIGHT_DT))*25/2, c="C5")
vbar(RIGHT_E, f"{RIGHT_E} keV", err=UC_E30, c="C6")
plt.xlabel("E (keV)")
plt.ylabel("Counts / hour")
plt.legend()
plt.show()


# histDT(sh28_dt3, label=f"Bare Stilbene (Unknown time ~ 28 hr), Ch3", scale=1, style="solid")
# sh284_counts, bins = histDT(sh28_dt4, label=f"37.5mm Fe CLYC (Unknown time ~ 28 hr), Ch4", scale=1, style="solid")

# histDT(bare_dt3, label="Bare Stilbene (72 hr), Ch 3", scale=1, style="dashed")
# bare4_counts, bins = histDT(bare_dt4, label="Bare CLYC (72 hr), Ch 4", scale=1, style="dashed")
# # vbar(LEFT_E, f"{LEFT_E} keV")
# # vbar(RIGHT_E, f"{RIGHT_E} keV")
# # plt.xlabel("E (keV)")
# plt.xlabel("dt (ns)")
# plt.ylabel("Counts")
# plt.legend()
# plt.show()

# plt.stairs(sh284_counts / bare4_counts, bins, label=" 37.5mm Fe / Bare (CLYC, Ch 4)")
# vbar(LEFT_E, f"{LEFT_E} keV")
# vbar(RIGHT_E, f"{RIGHT_E} keV")
# plt.xlabel("E (keV)")
# plt.ylabel("Ratio of count rates")
# plt.legend()
# plt.show()


def load_mockdt(path):
    '''T spectra are loaded in histogram/spectrum form already,
    so to rebin in coarser bins (this doesn't work for finer bins),
    we simply add one 'value' of halfway between the two bins to a new list of mock data'''
    # path = f"dat_Bare/raw2_ch{channel_num}_TOFSPEC.txt"
    with open(path, "r") as tofspec:
        text = tofspec.read()
        counts = [int(s) for s in text.split("\n") if not len(s)==0]
        bins = np.linspace(-1, 1, len(counts)+1) # start with just a scale of 0 to 1 then rescale after building a new series of mock data
        db = bins[1] - bins[0]
        vals = [] # mock data
        for i in range(len(counts)):
            vals.extend([bins[i]+0.5*db] * counts[i])
        
        return scale_tofspec(np.array(vals))

bareraw_dt4 = load_mockdt("dat_Bare/raw2_ch4_TOFSPEC.txt")
bareuf_dt4 = load_mockdt("dat_Bare/raw_ch4_TOFSPEC.txt")

sh22raw_dt4 = load_mockdt("dat_22hr/raw_ch4_TOFSPEC.txt")
sh22raw_dt0 = load_mockdt("dat_22hr/raw_ch0_TOFSPEC.txt")

MIN_DT = -30*2
MAX_DT = 130*2
NUM_BINS = 100
histDT(sh22_dt4, "Shielded CLYC, 72 hr")
# histDT(bareuf_dt4, "Unfiltered Bare CLYC, 72 h")
# histDT(sh22raw_dt0, "Raw Shielded")
# histDT(sh22raw_dt4, "")
plt.xlabel("dt (ns)")
plt.ylabel("Counts")
plt.legend()
plt.show()

