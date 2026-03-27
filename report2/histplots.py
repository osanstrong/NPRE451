import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm, poisson, binom, linregress
from math import floor, ceil

dat_20x30s:list = [
383, 
379,
340,
351,
374,
366,
369,
350,
374,
380,
357,
359,
364,
364,
357,
341,
369,
359,
345,
411,
]
# print(dat_20x30s)

dat_200x5s:list = [
55, 
52,
63,
64,
46,
59,
54,
59,
48,
54,
65,
68,
61,
59,
68,
53,
56,
51,
51,
56,
48,
52,
64,
53,
66,
65,
53,
71,
63,
57,
63,
55,
56,
55,
49,
58,
58,
55,
59,
66,
64,
51,
67,
58,
58,
70,
72,
46,
71,
61,
47,
62,
57,
51,
66,
59,
55,
54,
64,
52,
57,
57,
58,
52,
66,
45,
56,
65,
71,
60,
51,
63,
60,
52,
66,
49,
53,
58,
62,
58,
55,
75,
55,
63,
61,
62,
62,
42,
55,
48,
63,
48,
66,
70,
66,
53,
58,
62,
51,
73,
68,
48,
67,
47,
63,
59,
60,
60,
58,
50,
58,
56,
66,
64,
56,
58,
63,
75,
58,
41,
57,
63,
57,
58,
67,
55,
55,
57,
44,
67,
65,
66,
55,
44,
55,
73,
58,
61,
66,
61,
59,
63,
64,
58,
53,
59,
60,
49,
68,
65,
53,
47,
67,
64,
59,
70,
57,
56,
59,
62,
61,
48,
59,
63,
61,
59,
65,
49,
44,
67,
63,
53,
52,
55,
56,
69,
53,
48,
68,
60,
66,
60,
54,
54,
65,
67,
56,
73,
48,
63,
69,
63,
58,
62,
52,
57,
59,
60,
46,
65,
]
# print(dat_200x5s)

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

# print(mean(dat_20x30s))

select_dat = dat_20x30s
dat_mean = mean(select_dat)
dat_stdev = sample_stddev(select_dat)

print(f"Data mean: {dat_mean}, variance: {dat_stdev**2}")

xmin = dat_mean - dat_stdev*4
xmax = dat_mean + dat_stdev*4

# Now make bins manually
nb = 20
bins = np.linspace(xmin, xmax, nb)

eval_per_bin = 100
more_bins = np.linspace(xmin, xmax, nb*eval_per_bin)

bw = (xmax-xmin)/nb
norm_scl = len(select_dat)*bw



def scl_to_data(fit_vals):
    # normalize so distrib sums to 1
    fit_sum = sum(fit_vals) * bw
    print(fit_sum)
    new_fit_vals = fit_vals / fit_sum
    
    # fit to new sum
    new_fit_vals *= len(select_dat)
    print(sum(new_fit_vals) * bw)
    return new_fit_vals * bw



counts, bins = np.histogram(select_dat, bins=bins)
print(f"{len(counts)} vs {len(bins)}")
x = np.linspace(xmin, xmax, 250)
# k = range(ceil(xmin), floor(xmax))

k = [int((more_bins[i]+more_bins[i+1])*0.5) for i in range(len(more_bins)-1)] #halfway through each bin, rough approximation
k = [int(b) for b in more_bins[1:]]
k = np.arange(floor(xmin), ceil(xmax))

def rediscretize(x_vals, counts, new_bins):
    new_counts = np.zeros(len(new_bins)-1)

    j = 0 # index of current bin's left side
    #For each x
    for i in range(len(x_vals)):
        x = x_vals[i]
        c = counts[i]
        b_left = new_bins[j]
        b_right = new_bins[j+1]

        x_left = x-0.5
        x_right = x+0.5

        p_contained = min(x_right, b_right) - max(x_left, b_left)
        new_counts[j] += p_contained * c
        
        # If touching the next bin, also move on 
        if x_right > b_right:
            j += 1
            if j >= len(new_counts):
                break
            b_left = new_bins[j]
            b_right = new_bins[j+1]
            
            p_contained = min(x_right, b_right) - max(x_left, b_left)
            new_counts[j] += p_contained * c
    
    return new_counts



        
print(k)

# binom fit is more annoying
# μ = np
# s2 = np(1-p) = μ(1-p)
# 0 = np - np² - s2
# n = μ/p
# 0 = μ - μp - s2
# μp= μ - s2
# p = (μ-s2) / μ
# p = (dat_mean - dat_stdev**2) / dat_mean
# n = int(dat_mean/p)
n = int(dat_mean**2 / (dat_mean - dat_stdev**2))
p = dat_mean / n #Fix to fit the mean while keeping integer 


# fit_norm = norm.pdf(x,dat_mean, dat_stdev) * norm_scl
fit_norm = scl_to_data(norm.pdf(k, dat_mean, dat_stdev))
print(f"K used: {k}")
fit_pois = scl_to_data(poisson.pmf(k, dat_mean))
fit_bino = scl_to_data(binom.pmf(k, n, p))

sum_norm = rediscretize(k, fit_norm, bins)
print(sum_norm)
sum_pois = rediscretize(k, fit_pois, bins)
print(fit_pois)
sum_bino = rediscretize(k, fit_bino, bins)
# sum_pois, bins = np.histogram(fit_pois, bins=bins)
# sum_bino, bins = np.histogram(fit_bino, bins=bins)
# print(f"fit binom: {dat_mean/p}, {p}")

# fit_norm = scl(fitdw_norm, n**2)
# fit_pois = scl(fit_pois, n**2)
# fit_bino = scl(fit_bino, n**2)

# mu, std = norm.fit(select_dat)
# print(f"Manual mean std: {dat_mean}+- {dat_stdev}")
# print(f"Vs scipy of {mu}+- {std}")

distributions:dict = {
    "Experimental Trials": counts,
    # "Poisson Fit": fit_pois,
    # "Normal Fit": fit_norm,
    # "Binomial Fit": fit_bino,
    "Poisson Fit": sum_pois,
    "Normal Fit": sum_norm,
    "Binomial Fit": sum_bino,
}



# plt.stairs(counts, bins, label="Actual data")
# plt.plot(x, fit_norm, label="Normal fit")
# plt.stairs(fit_pois, k, label="Poisson Fit")
# plt.stairs(fit_bino, k, label="Binomial Fit")

# plt.step(bins[1:], [c for c in counts], label="Actual data")
# plt.step(bins[1:], fit_norm, label="Normal Fit")
# plt.step(bins[1:], fit_pois, label="Poisson Fit")
# plt.step(bins[1:], fit_bino, label="Binomial Fit")
# plt.hist(bins[:-1], bins, weights=fit_bino, fill=False, label="Binomial fit")
# plt.hist(bins[:-1], bins, weights=fit_pois, fill=False, label="Poisson fit")
# plt.hist(bins[:-1], bins, density=True, weights=counts, fill=False, label="Actual data")


for d in distributions:
    plt.fill_between(bins[1:], distributions[d], step="pre", alpha=0.3)

for d in distributions:
    plt.step(bins[1:], distributions[d], label=d)   


plt.xlabel("Counts in Trial")
plt.ylabel("Number of Trials")
plt.legend()
plt.show()


#### Plot Manual entry of a spectrogram plot

# manual_bins = [794.5, 795.5, 796.5, 797.5]
# vals = [41, 3537, 19]
# plt.hist(manual_bins[:-1], manual_bins, weights=vals, label="Acquired Spectrum at 4V")
# plt.xlabel("Bin (Marker #)")
# plt.ylabel("Counts")
# plt.legend()
# plt.show()


#### Plot V vs Spectrum Bin

v_vals = [4.0, 5.063, 6.000, 7.188]
c_vals = [759.96, 1024.03, 1207.32, 1441.54]

regr = linregress(v_vals, c_vals)

def plot_reg(res):
    xmin, xmax = plt.xlim()
    y = [xmin*res.slope + res.intercept, xmax*res.slope + res.intercept]
    plt.plot([xmin, xmax], y, 'r', label=f"Linear Fit: Slope {res.slope:.3f}, Intercept {res.intercept:.3f}, R^2:{res.rvalue**2:.3f}")
    # print("Fit")
    plt.xlim(xmin, xmax)

# plt.scatter(v_vals, c_vals, label="Recorded Spectra")
# plot_reg(regr)
# plt.xlabel("Voltage (V)")
# plt.ylabel("Centroid Peak (Channel #)")
# plt.legend()
# plt.show()