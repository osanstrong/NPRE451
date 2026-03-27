import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stat


# Experiment 1
biases = [
600, 
610,
620,
630,
640,
650,
660,
670,
680,
690,
700,
710,
720,
730,
740,
750,
800,
850,
900,
950,
1000
]

counts = [
00 , 
00 , 
00 , 
10 ,
37 ,
65 ,
244, 
372,
464,
493,
610,
523,
588,
612,
603,
611,
681,
693,
700,
740,
1000
]

plat_idx = [10, -2]
plat_v = biases[10:-1]
plat_c = counts[10:-1]

reg = stat.linregress(plat_v, plat_c)
m = reg.slope
b = reg.intercept

plat_creg = np.array(plat_v)*m + b

slope = (counts[plat_idx[-1]]-counts[-plat_idx[0]]) / (biases[plat_idx[-1]] - biases[plat_idx[0]])


# plt.xlabel("High voltage power supply Bias setting (V)")
# plt.ylabel("Measured Counts over 30 seconds")
# plt.plot(biases, counts, label="Measured Counts")
# # plt.plot([biases[i] for i in plat_idx], [counts[i] for i in plat_idx], label=f"Plateau Region (Slope: {slope} Counts/Volt)")
# plt.plot(plat_v, plat_creg, label=f"Plateau region (Slope: {m:.3f} counts per V)")

# plt.legend()
# plt.show()


#### Experiment 3
trials = [i for i in range(1, 34)]
counts = [
260, 
305,
298,
297,
295,
304,
261,
305,
293,
264,
302,
286,
249,
304,
244,
232,
260,
256,
246,
241,
245,
243,
231,
224,
231,
214,
192,
201,
201,
208,
176,
164,
178,
]

times = [90*i + 60 for i in range(0, 33)] #Approximate time elapsed in seconds

counts = np.array(counts)
times = np.array(times) / 60 #Convert to minutes

res = stat.linregress(times, np.log(counts))
print(f"Linear regression: {res}")

m = res.slope
b = res.intercept
fit_counts = np.exp(times*m + b)

hl = -np.log(2) / m
print(f"Half life: {hl} minutes")

plt.plot(times, counts, label="Measured counts of irradiated indium")
plt.plot(times, fit_counts, label="Exponential fit")
bot, top = plt.ylim()
plt.ylim(0, top)

plt.xlabel("Time elapsed (m)")
plt.ylabel("Counts detected in 60-second interval")

plt.legend()
plt.show()

