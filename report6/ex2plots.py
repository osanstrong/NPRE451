import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
# import scipy.stats as stat
import statsmodels.api as sm

al_x = np.array([25.4*i for i in range(5)])
al_r = [
    4.81,
    3.23,
    2.34,
    1.75,
    0.92,
]

fe_x = np.array([12.5*i for i in range(5)])
fe_r = [
    4.81,
    2.55,
    1.39,
    0.57,
    0.23,
]

pb_x = np.array([7.2*i for i in range(5)])
pb_r = [
    4.81,
    2.66,
    2.03,
    1.18,
    0.87
]

base_E = 1332.5

al_C = [
    1443, 74,
    967, 63,
    703, 54,
    525, 47,
    275, 40
]
fe_C = [
    1443, 74,
    763, 59,
    418, 47,
    170, 40,
    68, 32,
]
pb_C = [
    1443, 74,
    797, 63,
    608, 51,
    353, 45,
    260, 37,
]

al_y = (
    np.array([al_C[i*2] for i in range(5)]),
    np.array([al_C[i*2]-al_C[i*2+1] for i in range(5)]),
    np.array([al_C[i*2]+al_C[i*2+1] for i in range(5)]),
    np.array([al_C[i*2+1]**2 for i in range(5)])
)
al_yl = [np.log(y) for y in al_y]

# res = stat.linregress(al_x, al_yl[0])
al_xC = sm.add_constant(al_x)
print(f"al_x: {al_x}")
model = sm.WLS(np.log(al_y[0]), al_xC, weights=al_y[3])
res = model.fit()
coeffs = res.params
I0 = np.exp(coeffs[0])
mu = -coeffs[1]
print(f"Fit: I₀={I0:.0f} counts, μ={mu:.5f}mm⁻¹")
plt.scatter(al_x, al_y[0], label="Attenuated by aluminum")
plt.vlines(al_x, al_y[1], al_y[2], label="Uncertainty")
plt.plot(al_x, I0*np.exp(-mu*al_x), label=f"WLS Fit: I₀={I0:.0f} counts, μ={mu:.5f}mm⁻¹", c="C1")
plt.yscale("log")
plt.xlabel("Distance in material (mm)")
plt.ylabel("ROI Net Area (Counts)")
plt.legend()
plt.show()


fe_y = (
    np.array([fe_C[i*2] for i in range(5)]),
    np.array([fe_C[i*2]-fe_C[i*2+1] for i in range(5)]),
    np.array([fe_C[i*2]+fe_C[i*2+1] for i in range(5)]),
    np.array([fe_C[i*2+1]**2 for i in range(5)])
)
fe_yl = [np.log(y) for y in fe_y]

fe_xC = sm.add_constant(fe_x)
print(f"fe_x: {fe_x}")
model = sm.WLS(np.log(fe_y[0]), fe_xC, weights=fe_y[3])
res = model.fit()
coeffs = res.params
I0 = np.exp(coeffs[0])
mu = -coeffs[1]
print(f"Fit: I₀={I0:.0f} counts, μ={mu:.5f}mm⁻¹")
plt.scatter(fe_x, fe_y[0], label="Attenuated by iron")
plt.vlines(fe_x, fe_y[1], fe_y[2], label="Uncertainty")
plt.plot(fe_x, I0*np.exp(-mu*fe_x), label=f"WLS Fit: I₀={I0:.0f} counts, μ={mu:.5f}mm⁻¹", c="C1")
plt.yscale("log")
plt.xlabel("Distance in material (mm)")
plt.ylabel("ROI Net Area (Counts)")
plt.legend()
plt.show()



pb_y = (
    np.array([pb_C[i*2] for i in range(5)]),
    np.array([pb_C[i*2]-pb_C[i*2+1] for i in range(5)]),
    np.array([pb_C[i*2]+pb_C[i*2+1] for i in range(5)]),
    np.array([pb_C[i*2+1]**2 for i in range(5)])
)
pb_yl = [np.log(y) for y in pb_y]

pb_xC = sm.add_constant(pb_x)
print(f"pb_x: {pb_x}")
model = sm.WLS(np.log(pb_y[0]), pb_xC, weights=pb_y[3])
res = model.fit()
coeffs = res.params
I0 = np.exp(coeffs[0])
mu = -coeffs[1]
print(f"Fit: I₀={I0:.0f} counts, μ={mu:.5f}mm⁻¹")
plt.scatter(pb_x, pb_y[0], label="Attenuated by lead")
plt.vlines(pb_x, pb_y[1], pb_y[2], label="Uncertainty")
plt.plot(pb_x, I0*np.exp(-mu*pb_x), label=f"WLS Fit: I₀={I0:.0f} counts, μ={mu:.5f}mm⁻¹", c="C1")
plt.yscale("log")
plt.xlabel("Distance in material (mm)")
plt.ylabel("ROI Net Area (Counts)")
plt.legend()
plt.show()

al_mu = 0.1498 # inv cm
fe_mu = 0.5814 # 
pb_mu = 0.6055

al_rho = 2.699 # g / cm3
fe_rho = 7.874 # g / cm3
pb_rho = 11.35 # g / cm3

al_ma = al_mu / al_rho
fe_ma = fe_mu / fe_rho
pb_ma = pb_mu / pb_rho

print(f"Density (g/cm³) & Mass attenuation coefficients (cm²/g):")
print(f"al: {al_rho:.3f} {al_ma:.4f}")
print(f"fe: {fe_rho:.3f} {fe_ma:.4f}")
print(f"pb: {pb_rho:.3f} {pb_ma:.4f}")

# cm²/g @ 1.25 MeV and 1.5 MeV
al_NISTma = [0.05496, 0.05006]
fe_NISTma = [0.05350, 0.04883]
pb_NISTma = [0.05876, 0.05222]
def lerp(ab):
    a, b = ab
    t = (1.3325-1.25) / (1.5-1.25)
    return (1-t) * a + t * b
print("NIST Mass attenuation coeffs (cm²/g):")
print(f"al: {lerp(al_NISTma):.4f}")
print(f"fe: {lerp(fe_NISTma):.4f}")
print(f"pb: {lerp(pb_NISTma):.4f}")