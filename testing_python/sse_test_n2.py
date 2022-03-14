import numpy as np
import numpy.random as npr

from sse import *

N = 8
Nb = N
M = 12

T = 1.05
beta = 1 / T
rng = npr.default_rng(npr.MT19937(seed=5))

spin = np.array([1, 1, -1, -1, 1, -1, 1, -1])

bond_site = np.zeros((2, N), dtype=np.int64)
for b in range(N):
    bond_site[0, b] = b
    bond_site[1, b] = np.mod(b + 1, N) 
print(f"Bond Site = {bond_site[0, :]}")
print(f"            {bond_site[1, :]}")


opstring = np.array([4, 2, 9, 13, 6, 0, 0, 4, 13, 0, 9, 14], dtype=np.int64)
n = np.count_nonzero(opstring)
vertex_init, vertex_last, vertex_list = create_vertex_list(M, N, opstring, bond_site)

print(f"Spin State = {spin}")
print(f"n = {n}")
print("p | a[p]  |  b[p]  |  opstring[p] ")
for p in range(M):
    a = 0
    b = 0
    if opstring[p] != 0:
        a = np.mod(opstring[p], 2) + 1
        b = opstring[p] // 2 - 1
    print(f"{p} |   {a}   |   {b}    |    {opstring[p]}")
print(f"Vertex list = ")
for p in range(M):
    print(f"[{4 * p}] {vertex_list[4 * p]},  " + 
          f"[{4 * p + 1}] {vertex_list[4 * p + 1]},  " +
          f"[{4 * p + 2}] {vertex_list[4 * p + 2]},  " +
          f"[{4 * p + 3}] {vertex_list[4 * p + 3]}]")
print()

spin, opstring, n = diag_update(opstring, n, M, spin, bond_site, Nb, beta, rng)
print(f"Spin State = {spin}")
print(f"n = {n}")
print("p | a[p]  |  b[p]  |  opstring[p] ")
for p in range(M):
    a = 0
    b = 0
    if opstring[p] != 0:
        a = np.mod(opstring[p], 2) + 1
        b = opstring[p] // 2 - 1
    print(f"{p} |   {a}   |   {b}    |    {opstring[p]}")
print()


# opstring, spin = loop_update(M, N, opstring, spin, vertex_list, vertex_init, rng)
# print(f"Spin State = {spin}")
# print(f"n = {np.count_nonzero(opstring)}")
# print("p | a[p]  |  b[p]  |  opstring[p] ")
# for p in range(M):
#     a = 0
#     b = 0
#     if opstring[p] != 0:
#         a = np.mod(opstring[p], 2) + 1
#         b = opstring[p] // 2 - 1
#     print(f"{p} |   {a}   |   {b}    |    {opstring[p]}")
# print()



