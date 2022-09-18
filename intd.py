import numpy as np

f = np.loadtxt("data.txt", float)

N = len(f) - 1
a = f[0, 0]
b = f[N, 0]
h = (b - a) / N

s_Tr = (f[0, 1] + f[N, 1]) / 2
for i in range(1, N):
    s_Tr += f[i, 1]
I_Tr = h * s_Tr

s_Si = f[0, 1] + f[N, 1]
for i in range(1, N):
    if i % 2 == 0:
        s_Si += 2 * f[i, 1]
    else:
        s_Si += 4 * f[i, 1]
I_Si = h / 3 * s_Si

print("Trapezoidal Integral = ", I_Tr)
print("Simpson Integral = ", I_Si)
