import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

T = 180 # 180 seconds for every measurement

vac_dist = [1, 17, 25, 33, 49] #mm
vac_x = [895.1, 856.94, 815.68, 792.3, 666.6]
vac_I = [g/T for g in [
    2.824395e6, 
    8.399e5, # or 8.522e5
    3.405e5,
    2.709e5,
    1.272e5,
]]

E = 5486 #keV
C = E / vac_x[0] # E corresponding to each bin

A = 450 #mm²
PI = np.pi
a = 2*np.sqrt(A / PI) # mm
def sangle(d):
    return (2*PI*(1 - d/(np.hypot(d, a))))
# print([sangle(d) for d in vac_dist])
def norm(I, d):
    return 4*PI*I / sangle(d)
    # return I
norm_I = [4*PI*vac_I[i]/sangle(vac_dist[i]) for i in range(len(vac_dist))]

# plt.scatter(vac_dist, vac_I, label="In vacuum")
# plt.scatter(vac_dist, norm_I, label="In vacuum, normalized by solid angle")
# plt.xlabel("Distance (mm)")
# # plt.ylabel("Integrated peak intensity (counts)")
# plt.ylabel("Peak intensity (counts per second)")
# plt.legend()
# plt.show()


air_dist = [1, 17, 25, 33, 49] #mm
air_I = [g/T for g in [
    2.841e6,
    8.546e5,
    3.44e5,
    2.72e5,
    1
]]
air_x = [x*C for x in [
    865.58,
    661.38,
    438.20,
    331.65,
    0
]]
air_FWHM = [f*C for f in [
    22.22,
    60.39,
    55.56,
    58.84,
    0
]]

# plt.scatter(air_dist, air_I, label="In air")
# plt.scatter(air_dist, [norm(air_I[i], air_dist[i]) for i in range(5)], label="In air, normalized by solid angle")
# plt.xlabel("Distance (mm)")
# # plt.ylabel("Integrated peak intensity (counts)")
# plt.ylabel("Peak intensity (counts per second)")
# plt.legend()
# plt.show()

# plt.scatter(air_dist[:4], air_x, label="In air")
# plt.xlabel("Distance (mm)")
# plt.ylabel("Peak Energy (keV)")
# plt.legend()
# plt.show()
# plt.scatter(air_dist[:4], air_FWHM, label="In air")
# plt.xlabel("Distance (mm)")
# plt.ylabel("Peak FWHM (keV)")
# plt.legend()
# plt.show()

al_thickness = 18 #microns
al_dist = [al_thickness * n for n in range(1,6)]
al_dat = [
    [610.5, 60.19, 3195, 2542], # I in rates this time
    [293.68, 116.28, 2146, 1112],
    [8.51, 2.61, 58, 7],
    [8.77, 1.85, 17, 10],
    [24, 1, 0.15, 0]
]
al_x = [d[0]*C for d in al_dat]
al_f = [d[1]*C for d in al_dat]
al_i = [d[2] for d in al_dat] # For now using gross

# plt.scatter(al_dist, al_x, label="Aluminum")
# plt.xlabel("Thickness of material (microns)")
# plt.ylabel("Peak location (keV)")
# plt.legend()
# plt.show()

# plt.scatter(al_dist, al_f, label="Aluminum")
# plt.xlabel("Thickness of material (microns)")
# plt.ylabel("Peak FWHM (keV)")
# plt.legend()
# plt.show()

# plt.scatter(al_dist, al_i, label="Aluminum")
# plt.xlabel("Thickness of material (microns)")
# plt.ylabel("Peak intensity (counts per second)")
# plt.legend()
# plt.show()

#### Mylar ####
my_thickness = 3.6 #microns
my_dist = [my_thickness * n for n in range(1,6)]
my_dat = [
    [758.58, 26.58, 3328, 2973],
    [673, 38.62, 3308, 2701],
    [491, 46.84, 1114, 833], # 3 left peak
    # [589.7, 55.46, 441, 330], # 3 right peak
    [512.82, 57.72, 1791, 1199],
    [434.42, 65.60, 2715, 1771]
]
my_x = [d[0] for d in my_dat]
my_f = [d[1]*C for d in my_dat]
my_i = [d[2]*C for d in my_dat] # For now using gross

# plt.scatter(my_dist, my_x, label="Mylar")
# plt.xlabel("Thickness of material (microns)")
# plt.ylabel("Peak location (keV)")
# plt.legend()
# plt.show()

# plt.scatter(my_dist, my_f, label="Mylar")
# plt.xlabel("Thickness of material (microns)")
# plt.ylabel("Peak FWHM (keV)")
# plt.legend()
# plt.show()

# plt.scatter(my_dist, my_i, label="Mylar")
# plt.xlabel("Thickness of material (microns)")
# plt.ylabel("Peak intensity (counts per second)")
# plt.legend()
# plt.show()


#### Get peak energy of simulated transmissions

BW = 1e4 #eV
min_E = 0.1e7
max_E = 0.7e7
NB = int(max_E / BW)
BINS = np.linspace(min_E, max_E, NB)

N_P = 1000


# Import and test
def get_spectrum(fpath):
    # with open(fpath, 'r') as file:
    #     print(file.read())
    df = pd.read_csv(fpath, sep="\s+", skiprows=11)
    vals = df['(eV)']
    vals, bins = np.histogram(vals, BINS)
    return vals*1000
    # return df['(eV)']
    # return df
def get_peak_E_hist(c, bins):
    return bins[np.argmax(c)]
def get_peak_E(fpath):  
    c= get_spectrum(fpath)
    return get_peak_E_hist(c, BINS)
# c, bins = get_spectrum("data/Al, 3, TRANSMIT.txt")
# print(get_peak_E_hist(c, bins))
# plt.step(BINS[1:], c, label="Simulated 3 layer of Al")
# plt.xlabel("E (eV)")
# plt.ylabel("Number of counts")
# plt.legend()
# plt.show()
# print(get_spectrum("data/Al, 1, TRANSMIT.txt")

al_Sp = [get_spectrum(path)/1e3 for path in [f"data/Al, {n}, TRANSMIT.txt" for n in range(1,6)]]
my_Sp = [get_spectrum(path)/1e3 for path in [f"data/Mylar, {n}, TRANSMIT.txt" for n in range(1,6)]]
air_TRIMdist = [10, 20, 30]
air_Sp = [get_spectrum(path)/1e3 for path in [f"data/Air, {d}, TRANSMIT.txt" for d in air_TRIMdist]]

# for i in range(5):
#     plt.step(BINS[1:]/1e6, al_Sp[i], label=f"Aluminum, {al_dist[i]} microns")
# for i in range(5):
#     plt.step(BINS[1:]/1e6, my_Sp[i], label=f"Mylar, {my_dist[i]} microns")
# for i in range(3):
#     plt.step(BINS[1:]/1e6, air_Sp[i], label=f"Aluminum, {air_TRIMdist[i]} mm")
# plt.xlabel("E (MeV)")
# plt.ylabel("Number of counts")
# plt.legend()
# plt.show()


al_E = [get_peak_E(path)/1e3 for path in [f"data/Al, {n}, TRANSMIT.txt" for n in range(1,6)]]
# al_E[2] = (al_E[1] + al_E[3])/2
# plt.scatter(al_dist, al_E, label="Aluminum")
# plt.xlabel("Thickness of material (microns)")
# plt.ylabel("Peak energy (keV)")
# plt.legend()
# plt.show()

my_E = [get_peak_E(path)/1e3 for path in [f"data/Mylar, {n}, TRANSMIT.txt" for n in range(1,6)]]
# plt.scatter(my_dist, my_E, label="Mylar")
# plt.xlabel("Thickness of material (microns)")
# plt.ylabel("Peak energy (keV)")
# plt.legend()
# plt.show()
air_E = [get_peak_E(path)/1e3 for path in [f"data/Air, {d}, TRANSMIT.txt" for d in [10, 20, 30]]]



# my_E = my_x
my_dEdx = [(my_E[i] - my_E[i+1])/(my_dist[i+1]-my_dist[i]) for i in range(len(my_dist)-1)]
my_Ed = [(my_E[i] + my_E[i+1])/2 for i in range(len(my_E)-1)]
# plt.scatter(my_Ed, my_dEdx, label="Mylar")
# print(my_dEdx)
# plt.xlabel("Energy (keV)")
# plt.ylabel("Stopping power -dE/dx (keV/micron)")
# plt.legend()
# plt.show()

# air_E = air_x


air_E = [get_peak_E(f"data/Air, {d}, TRANSMIT.txt")/1e3 for d in air_dist[:-1]] + [0]

air_dEdx = [(air_E[i] - air_E[i+1])/(air_dist[i+1]-air_dist[i]) for i in range(len(air_dist)-1)]
air_Ed = [(air_E[i] + air_E[i+1])/2 for i in range(len(air_E)-1)]
# plt.scatter(air_Ed, air_dEdx, label="Air")
# print(air_dEdx)
# plt.xlabel("Energy (keV)")
# plt.ylabel("Stopping power -dE/dx (keV/mm)")
# plt.legend()
# plt.show()

# al_E = al_x

# al_dEdx = [(al_E[i] - al_E[i+1])/(al_dist[i+1]-al_dist[i]) for i in range(len(al_dist)-1)]
# al_Ed = [(al_E[i] + al_E[i+1])/2 for i in range(len(al_E)-1)]
# plt.scatter(al_Ed, al_dEdx, label="Aluminum")
# print(al_dEdx)
# plt.xlabel("Energy (keV)")
# plt.ylabel("Stopping power -dE/dx (keV/micron)")
# plt.legend()
# plt.show()


def load_dEdx(path, skip):
    df = pd.read_csv(path, delimiter="\s+", skiprows=skip)
    cols = df.columns
    foot = 13
    E = df[cols[0]][:-foot]
    E = E.astype(float)
    E[-27:] = E[-27:]*1000
    E /= 1000
    Se = df[cols[2]][:-foot].astype(float)
    Sn = df[cols[3]][:-foot].astype(float)
    St = Se + Sn
    print(E)
    return E, Se, Sn, St

def plot_dEdx(path, skip, mat, trans_est):
    # df = pd.read_csv(path, delimiter="\s+", skiprows=24)
    # cols = df.columns
    # E = df[cols[0]]
    # Se = df[cols[2]]
    # Sn = df[cols[3]]
    # St = Se + Sn
    E, Se, Sn, St = load_dEdx(path, skip)

    plt.plot(E, Se, label=f"Electric, in {mat}")
    plt.plot(E, Sn, label=f"Nuclear, in {mat}")
    plt.plot(E, St, label=f"Total, in {mat}")
    if not trans_est is None:
        plt.scatter(trans_est[0], trans_est[1], label=f"Total in {mat}, Transmission estimate")
    plt.xlabel("E (MeV)")
    plt.ylabel("Stopping Power -dE/dx (keV/micron)")
    plt.loglog()
    plt.legend()
    plt.show()

    return E, Se, Sn, St

plot_dEdx("data/SRIM_ Helium in Aluminum", 24, "aluminum", None)
plot_dEdx("data/SRIM_ Helium in Mylar", 26, "Mylar", ([e/1e3 for e in my_Ed], [e for e in my_dEdx]))
plot_dEdx("data/SRIM_ Helium in Air, Dry (gas)", 28, "air", ([e/1e3 for e in air_Ed], [e/1e3 for e in air_dEdx]))

def load_range(path, skip):
    df = pd.read_csv(path, delimiter="\s+", skiprows=skip)
    cols = df.columns
    foot = 13
    E = df[cols[0]][:-foot]
    E = E.astype(float)
    E[-27:] = E[-27:]*1000
    E /= 1000
    range = df[cols[4]][:-foot].astype(float)
    scale = df[cols[5]][:-foot].replace("um", 1).replace("mm", 1000).replace("A", 1e-4).astype(float)
    return E, range*scale

# air_r = load_range("data/SRIM_ Helium in Air, Dry (gas)", 28)
# al_r = load_range("data/SRIM_ Helium in Aluminum", 24)
# my_r = load_range("data/SRIM_ Helium in Mylar", 26)

# plt.plot(air_r[0], air_r[1]/1000, label="Air, SRIM")
# plt.errorbar([5.486], [21], yerr=[4], label="Air, observed")
# plt.xlabel("E (MeV)")
# plt.ylabel("Range (mm)")
# plt.legend()
# plt.show()

# plt.plot(al_r[0], al_r[1], label="Aluminum, SRIM")
# plt.errorbar([5.486], [36], yerr=[9], label="Aluminum, observed")
# plt.plot(my_r[0], my_r[1], label="Mylar, SRIM")
# plt.errorbar([5.486], [9], yerr=[1.8], label="Mylar, observed")
# plt.xlabel("E (MeV)")
# plt.ylabel("Range (microns)")
# plt.legend()
# plt.show()