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
    
    



print(f"Lengths: 0: {len(ch0)}, 3: {len(ch3)}, 4: {len(ch4)}")



def hist(dat, label):
    bins = np.linspace(min(dat), max(dat), 1000)
    counts, bins = np.histogram(dat, bins)
    plt.stairs(counts, bins, label=label)
# hist(ch4, label="Channel 4")
hist(ch4_dt, label="Channel 4")
# hist(ch3_dt, label="Channel 3")

plt.xlabel("dt (ns)")
plt.ylabel("Counts")
plt.legend()
plt.show()