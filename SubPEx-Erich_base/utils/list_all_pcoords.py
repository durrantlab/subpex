import json
import os
import math
import glob
import os

seg_dirs = glob.glob("traj_segs/*")
seg_dirs.sort()

if os.path.exists("traj_seg.data.json"):
    data = json.load(open("traj_seg.data.json", "r"))
else:
    data = {}

amin = 1e10
amax = 0
for seg_dir in seg_dirs:
    pcoors = glob.glob(seg_dir + "/*/pcoord.txt")

    vals = []
    for pcoor in pcoors:
        if pcoor in data:
            vals.append(data[pcoor])
        else:
            lines = open(pcoor).readlines()
            if len(lines) > 1:
                val = lines[-1].split()[-1]
                val = math.floor(100 * float(val)) / 100
                vals.append(val)
                data[pcoor] = val

    vals.sort()
    if len(vals) > 0:
        if vals[0] < amin:
            amin = vals[0]
        if vals[-1] > amax:
            amax = vals[-1]
        vals = [str(v).ljust(4,"0") for v in vals]
        print seg_dir, " " + vals[0] + "-" + vals[-1] + " ", " ".join(vals)

json.dump(data, open("traj_seg.data.json", "w"))

print "Min:", amin, "; Max:", amax
