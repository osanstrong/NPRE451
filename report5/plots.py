import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def readlist(path):
    with open(path, 'r') as file:
        lines = file.read().splitlines()
        return np.array([float(n) for n in lines])
    
def vbar(x, label, c=None, style="dashed"):
    [ymin, ymax] = plt.ylim()
    plt.vlines(x, ymin, ymax, label=label, colors=c, linestyles=style)
    plt.ylim(ymin, ymax)
    

calpath_inorg = "dat/extra/cali_Inorg (ch1)/RAW/CH1@DT5730S_2263_EspectrumR_cali_Inorg_20260227_134813.txt"
calpath_org = "dat/extra/cali_Org (ch0)/RAW/CH0@DT5730S_2263_EspectrumR_cali_Org_4_20260227_134200.txt"
spec_00_in = readlist(calpath_inorg)
spec_00_og = readlist(calpath_org)
chans_in = np.array([i for i in range(len(spec_00_in))])
chans_og = np.array([i for i in range(len(spec_00_og))])

# plt.plot(chans_in, spec_00_in, label="0 degrees C1")
# plt.plot(chans_og, spec_00_og, label="calibration org")
# plt.legend()
# plt.show()

CS137_PHOTOPEAK = 661.657 # keV
# Calculate the peak compton? I wonder if organic just doesn't have any photopeak
# m_e = 0 #kg
# c = 3e8 #m/s
mec2 = 511 # keV, rest mass in units of energy
J_per_MeV = 1.602e-13
J_per_keV = J_per_MeV * 1e-3
keV_per_J = 1/J_per_keV


def get_compton_Ee(Eg, theta):
    # Eg in keV
    Ee = Eg - get_compton_Egp(Eg, theta)
    return Ee
def get_compton_Egp(Eg, theta):
    # Eg in keV
    Egp = Eg / (1 + (Eg/mec2)*(1-np.cos(theta)))
    return Egp
def get_theta(Eg, phi):
    alpha = Eg/mec2
    # TODO: Figure this out with internet or pencil and paper if time
# print(f"Peak compton from {CS137_PHOTOPEAK} keV photon: {get_compton_Ee(CS137_PHOTOPEAK, np.pi)} keV")
CS137_COMPEAK = get_compton_Ee(CS137_PHOTOPEAK, np.pi)


C_INORG = CS137_PHOTOPEAK / 2626 # Known photopeak
C_ORG = CS137_COMPEAK / 248 # Known Compton edge
# print(f"Calibration:\n{C_INORG} keV/chan for inorganic\n{C_ORG} keV/chan for organic")

E_in = chans_in * C_INORG
E_og = chans_og * C_ORG

# plt.plot(E_in, spec_00_in, label="Inorganic, calibration")
# vbar(CS137_PHOTOPEAK, f"Theoretical photopeak: {CS137_PHOTOPEAK} keV", c="C2")

# plt.plot(E_og, spec_00_og, label="Organic, calibration")
# vbar(CS137_COMPEAK, f"Theoretical Compton edge: {CS137_COMPEAK:.3f} keV", c="C3")
# plt.xlabel("E (keV)")
# plt.ylabel("Counts")
# plt.xlim(-10, 870)
# plt.legend()
# plt.show()


angles = ["00", "30", "60", "90"]
angrads = [float(ang)*np.pi/180 for ang in angles]
paths_in = [f"dat/{ang}_inorg.txt" for ang in angles]
paths_og = [f"dat/{ang}_org.txt" for ang in angles]

specs_in = [readlist(path) for path in paths_in]
specs_og = [readlist(path) for path in paths_og]

ideal_Ee = [get_compton_Ee(CS137_PHOTOPEAK, rad) for rad in angrads]
ideal_Egp = [get_compton_Egp(CS137_PHOTOPEAK, rad) for rad in angrads]
print(f"Ee: {ideal_Ee}")
print(f"Egp: {ideal_Egp}")
# for i in range(len(angles)):
#     ang = angles[i]
#     rad = angrads[i]
#     plt.plot(E_in, specs_in[i], label=f"Inorganic, {ang} degrees", alpha=0.5)
# for i in range(len(angles)):
#     ang = angles[i]
#     vbar(ideal_Ee[i], f"Predicted recoil electron E, {ang} degrees", c=f"C{i}", style="dashed")
# for i in range(len(angles)):
#     ang = angles[i]
#     vbar(ideal_Egp[i], f"Predicted scattered gamma E, {ang} degrees", c=f"C{i}", style="solid")
    
# plt.xlabel("E (keV)")
# plt.ylabel("Counts")
# plt.xlim(-10, 1070)

# plt.legend()
# plt.show()


# for i in range(len(angles)):
#     ang = angles[i]
#     rad = angrads[i]
#     # ideal_CE = get_compton_Ee(CS137_PHOTOPEAK, rad)
#     plt.plot(E_og, specs_og[i], label=f"Organic, {ang} degrees")
#     # vbar(ideal_CE, f"Predicted Compton Edge, {ang} degrees", c=f"C{i}")
# for i in range(len(angles)):
#     ang = angles[i]
#     vbar(ideal_Ee[i], f"Predicted recoil electron E, {ang} degrees", c=f"C{i}", style="dashed")
# for i in range(len(angles)):
#     ang = angles[i]
#     vbar(ideal_Egp[i], f"Predicted scattered gamma E, {ang} degrees", c=f"C{i}", style="dotted")
  
# plt.xlabel("E (keV)")
# plt.ylabel("Counts")
# plt.xlim(-10, 870)
# plt.legend()
# plt.show()