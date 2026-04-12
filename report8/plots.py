import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress as linreg

layers = np.array([i for i in range(4)])
l_conc = 10 # cm
l_conc_unc = 0.1
t_conc = layers * l_conc
t_conc_unc = layers * l_conc_unc

l_steel = 2 # cm
l_steel_unc = 0.1
t_steel = layers * l_steel
t_steel_unc = layers * l_steel_unc

l_pe = 5 # cm
l_pe_unc = 0.1
t_pe = layers * l_pe
t_pe_unc = layers * l_pe_unc

d = 87
a = 11
F = a**2 / d**2


D_conc = np.array([
    70.0, 41.6, 26.9, 17.0
]) # μrem
D_steel = np.array([
    70.0, 60.0, 50.9, 43.6
])
D_pe = np.array([
    70.0, 47.9, 29.2, 13.6
])
ln_conc = np.log(D_conc/F)
ln_steel = np.log(D_steel/F)
ln_pe = np.log(D_pe/F)

reg_conc = linreg(t_conc, ln_conc)
line_conc = reg_conc.slope * t_conc + reg_conc.intercept
XSr_conc = -reg_conc.slope
R_conc = reg_conc.rvalue**2

reg_steel = linreg(t_steel, ln_steel)
line_steel = reg_steel.slope * t_steel + reg_steel.intercept
XSr_steel = -reg_steel.slope
R_steel = reg_steel.rvalue**2

reg_pe = linreg(t_pe, ln_pe)
line_pe = reg_pe.slope * t_pe + reg_pe.intercept
XSr_pe = -reg_pe.slope
R_pe = reg_pe.rvalue**2
# print(f"linreg of steel: {reg_steel}")

# plt.scatter(t_conc, ln_conc, label="Concrete", marker="d", c="C0")
# plt.plot(t_conc, line_conc, '--', label="Concrete (linear regression)", c="C0")
# plt.scatter(t_steel, ln_steel, label="Steel", marker="*", c="C1")
# plt.plot(t_steel, line_steel, '--', label="Steel (linear regression)", c="C1")
# plt.scatter(t_pe, ln_pe, label="Polyethylene", marker="X", c="C2")
# plt.plot(t_pe, line_pe, '--', label="Polyethylene (linear regression)", c="C2")

def push_xsr_plot():
    plt.xlabel("Thickness of shielding (cm)")
    plt.ylabel("ln(D(t) * d²/a²)")
    plt.legend()
    plt.show()
# push_xsr_plot()

plt.scatter(t_conc, ln_conc, label="Concrete", marker="d", c="C0")
plt.plot(t_conc, line_conc, '--', label=f"Linear Regression, ΣR={XSr_conc:.3f} cm⁻¹, R²={R_conc:.3f}", c="C0")
push_xsr_plot()

plt.scatter(t_steel, ln_steel, label="Steel", marker="*", c="C1")
plt.plot(t_steel, line_steel, '--', label=f"Linear Regression, ΣR={XSr_steel:.3f} cm⁻¹, R²={R_steel:.3f}", c="C1")
push_xsr_plot()

plt.scatter(t_pe, ln_pe, label="Polyethylene", marker="X", c="C2")
plt.plot(t_pe, line_pe, '--', label=f"Linear Regression, XR={XSr_pe:.3f} cm⁻¹, R²={R_pe:.3f}", c="C2")
push_xsr_plot()

# No graphing these actually
D_conc_Z = np.array([
    78.6, 42.2, 28.6, 17.8
]) / D_conc
D_conc_YZ = np.array([
    89.9, 40.5, 26.0, 19.4
]) / D_conc
D_conc_Y = np.array([
    78.8, 41.2, 28.3, 19.5
]) / D_conc

print("XZ & " + " & ".join([f"{D:.3f}" for D in D_conc_Z]) + " \\\\")
print("XYZ & " + " & ".join([f"{D:.3f}" for D in D_conc_YZ]) + " \\\\")
print("XY & " + " & ".join([f"{D:.3f}" for D in D_conc_Y]) + " \\\\")

