import subprocess
import numpy as np

read_in = lambda S, h, delta: f"1, 6, {S}, {delta}, {h}, 0.25, PBC\n10000, 1000000, 20"

S = 0.5
delta = 1.0

h_vals_1 = np.arange(0.0, S*1.0, S*0.025)
h_vals_2 = np.arange(S*1.0, S*3.0, S*0.25)
h_vals_3 = np.arange(S*3.0, S*4.0, S*0.025)

h_vals = np.concatenate([h_vals_1, h_vals_2, h_vals_3, [S*4.0]])

for i, h in enumerate(h_vals):
    with open("read.in", "w") as file:
        file.write(read_in(S, np.round(h, 6), delta))

    subprocess.call(f"./prepare_job.sh -n 20 -r 5-0 -j L6_h{i}", shell=True)
