import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from scipy.stats import norm


PEAK_BIN = 38.72
PEAK_DT = 6.5
FWHM_BIN = 7.28
FWHM_DT = 2.2
scale_tofspec = lambda dt_b: PEAK_DT +  (dt_b-PEAK_BIN)*(FWHM_DT/FWHM_BIN) - 0.073

b5 = 32
b20 = 90
m = (20 - 5) / (b20 - b5)
b = 5 - m*b5

scale_tofspec = lambda dt_b: m*(dt_b-2.25) + b + 0.064

def load_mockdt(path):
    '''T spectra are loaded in histogram/spectrum form already,
    so to rebin in coarser bins (this doesn't work for finer bins),
    we simply add one 'value' of halfway between the two bins to a new list of mock data'''
    # path = f"dat_Bare/raw2_ch{channel_num}_TOFSPEC.txt"
    with open(path, "r") as tofspec:
        text = tofspec.read()
        counts = [int(s) for s in text.split("\n") if not len(s)==0]
        bins = np.array(range(len(counts)+1))
        return counts, scale_tofspec(bins)
        bins = np.linspace(-1, 1, len(counts)+1) # start with just a scale of 0 to 1 then rescale after building a new series of mock data
        db = bins[1] - bins[0]
        vals = [] # mock data
        for i in range(len(counts)):
            vals.extend([bins[i]+0.5*db] * counts[i])
        
        return scale_tofspec(np.array(vals)), bins
   

def fit_func(x,a,mu,sigma,c):
    """gaussian function used for the fit"""
    return a * norm.pdf(x,loc=mu,scale=sigma) + c

na22_dt = load_mockdt("na22_calib/Na22_DT_Manual.txt")
counts, bins = na22_dt
centers = bins[:-1] + 0.5*(bins[1]-bins[0])
p1,_ = curve_fit(fit_func, centers, counts)
x = np.array(centers)
y_fit = fit_func(x, *p1)
a, mu, sigma, c = p1

print(f"Mean: {mu}, sigma: {sigma}")
plt.stairs(counts, bins, label="$^{22}$Na positron calibration")
plt.plot(x, y_fit, label=f"Fit: {a:.2f}*Norm(\n\t$\mu$={mu:.3f} ns,\n\t$\sigma$={sigma:.3f} ns\n) + {c:.2f}")
plt.xlabel("dt (ns)")
plt.ylabel("Counts (-)")
plt.legend()
plt.show()




# def fit_func(x,a,mu,sigma,c):
#     """gaussian function used for the fit"""
#     return a * norm.pdf(x,loc=mu,scale=sigma) + c

#make up some normally distributed data and do a histogram
y = 2 * np.random.normal(loc=1,scale=2,size=1000) + 2
no_bins = 20
hist,left = np.histogram(y,bins=no_bins)
centers = left[:-1] + (left[1] - left[0])

#fit the histogram
p0 = [2,0,2,2] #starting values for the fit
p1,_ = curve_fit(fit_func,centers,hist,p0,maxfev=10000)

#plot the histogram and fit together
fig,ax = plt.subplots()
ax.hist(y,bins=no_bins)
x = np.linspace(left[0],left[-1],1000)
y_fit = fit_func(x, *p1)
ax.plot(x,y_fit,'r-')
# plt.show()