import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stat

path = "L1D/scope_2_1.csv"

dat = pd.read_csv(path)
# print(dat)

x = dat['x-axis'][1:].astype(float)
x_units = dat['x-axis'][0]
y = dat['1'][1:].astype(float)
y_units = dat['1'][0]


dat2 = pd.read_csv("L1D/scope_3_1.csv")
x2 = dat2['x-axis'][1:].astype(float)
x_units2 = dat2['x-axis'][0]
y2 = dat2['1'][1:].astype(float)
y_units2 = dat2['1'][0]


ids = [
    "5_2",
    "6_1",
    "6_2",
    "7_1",
    "7_2",
    "8_1",
    "8_2"
]

names = [
    "Only Preamplifier",
    "Unipolar, Lowest Gain",
    "(Preamplifier Output)",
    "Bipolar, Lowest Gain",
    "(Preamplifier Output)",
    "Unipolar, 200x Coarse Gain",
    "(Preamplifier Output)"
]

indices = [1, 2, 3, 4, 5, 6]

set_paths = [f"L1D/scope_{i}.csv" for i in ids]


fig, ax = plt.subplots()
plt.plot(x2, y2, label="Directly Acquired Pulse (Fall time: 1.14ms, Ampl.: 9.68V)")
plt.plot(x, y, label="Attenuated Pulse (Fall time: 1.13ms, Ampl.: 5.00V)")
# for i in indices:
#     id = ids[i]
#     name = names[i]
#     path_id = f"L1D/scope_{id}.csv"
#     channel = id.split("_")[1]

#     pdat = pd.read_csv(path_id)
#     x = pdat['x-axis'][1:].astype(float)
#     x_units = pdat['x-axis'][0]
#     y = pdat[channel][1:].astype(float)
#     y_units = pdat[channel][0]
#     plt.plot(x, y, label=name)
    
plt.xlabel(f"Time ({x_units}s)")
plt.ylabel(f"Measurement ({y_units}s)")



plt.title("")


plt.legend()
plt.show()


#### Pulse Amp vs Pulse Height setting

pulse_height = [
    1.74,
    2.6,
    3.8,
    5,
]

pulse_amp_amp_mV = [
    4.25,
    8.25,
    8.2,
    8.5,
]

pulse_pre_amp_mV = [
    2.925,
    5.5,
    5.875,
    5.925,
]

#### Pulse Amp vs Coarse Gain

coarse_gain = [
    10,
    20,
    50,
    100
]

coarse_amp_amp_mV = [
    4,
    7.825,
    17.3,
    32.625,
]
res_ph = stat.linregress(pulse_height, pulse_amp_amp_mV)
res_cg = stat.linregress(coarse_gain, coarse_amp_amp_mV)
def plot_reg(res, ax):
    xmin, xmax = plt.xlim()
    y = [xmin*res.slope + res.intercept, xmax*res.slope + res.intercept]
    ax.plot([xmin, xmax], y, 'r', label=f"Amplifier Fit: Slope {res.slope:.3f}, Intercept {res.intercept:.3f}, R^2:{res.rvalue**2:.3f}")
    # print("Fit")
    plt.xlim(xmin, xmax)

# fig, ax = plt.subplots()
# # plt.xlabel("Pulse Height Setting")
# plt.xlabel("Coarse Gain Setting")
# plt.ylabel("Pulse Amplitude (mV)")


# # ax.scatter(pulse_height, pulse_pre_amp_mV, label="Preamplifier Only")
# # ax.scatter(pulse_height, pulse_amp_amp_mV, label="Amplifier Output")
# # plot_reg(res_ph, ax)

# ax.scatter(coarse_gain, coarse_amp_amp_mV, label="Amplifier Output")
# plot_reg(res_cg, ax)

# plt.legend()
# plt.show()
