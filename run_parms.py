import subprocess
import numpy as np

read_in = lambda S, h, delta: f"1, 6, {S}, {delta}, {h}, 0.25, PBC\n100000, 1000000, 20"

S = 0.5
delta = 1.0

h_vals = np.arange(0.0, S*4.0 + S*0.025, S*0.025)

for i, h in enumerate(h_vals):
    with open("read.in", "w") as file:
        file.write(read_in(S, np.round(h, 6), delta))

    subprocess.call(f"./preprare_job.sh -n 20 -r 5-0 -j L6_h{i}", shell=True)
