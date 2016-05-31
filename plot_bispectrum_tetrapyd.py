import numpy as np
import math
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
xs1 = []
ys1 = []

filename1 = "LISW_tetrapyd.dat"

f1 = file(filename1)
data = np.loadtxt(f1)

