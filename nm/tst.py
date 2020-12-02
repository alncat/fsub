from mmtbx.nm import tools
from mmtbx_nm_ext import *
from cctbx.array_family import flex
from pprint import pprint
import time

if (__name__ == "__main__"):
    nmval = []
    nmval = tools.read_nmval_file()
    for i in range(len(nmval)):
        print i, nmval[i]
    x = flex.double(50*50, 0.0)
    nx = 0
    for i in range(6):
        for j in range(i+1):
            x[nx] = i+j
            nx += 1
    for i in range(6, 50):
        for j in range(6, i+1):
            x[nx] = i+j
            nx += 1
    t1 = time.time()
    s = unpack_x(x = x, n_modes = 50, zero_mode_flag = True)
    t2 = time.time()
    print t2 - t1
    pprint(list(s))
