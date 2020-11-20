import numpy as np
from scipy.integrate import ode
import os
#  set the 'MPLCONFIGDIR' to '/tmp' to get a writable environment for matplotlib
os.environ['MPLCONFIGDIR'] = '/tmp'
import matplotlib.pyplot as plt
import csv
from KINmodule import *

# time discretizzation
t0 = 0
tf = 30
dt = 0.05 # [s]
N   = round((tf-t0)/dt)+1
t = np.array(range(N))*dt+t0

# data loading
[u0, L, l, b] = dataBuilder()

# initializing matrix A: statedot system matrix (input dependent)
tmp = [[0.0 for i in range(7)] for j in range(7)]  # temporary array
A = np.array(tmp)
# constant part of matrix A
for i in range(6):
    A[(i + 1), 0] = b[i] / L
    A[(i + 1), (i + 1)] = -l[i]
# set the STATEDOT parameters (see KINmodule)
# REACTIVITY FUNCTION CHOICE HERE (first argument of param)
param = [rho_step, L, l, b,A]

######################################################################
# problem solution ###################################################

# SOLUTION initializzation
ts = []
sol = []
# SOLUTOR parameters (same solutor as MATLAB ode15s )
ode15s = ode(statedot)
# max step limitation necessary
ode15s.set_integrator('zvode', method='bdf', order=15, max_step=dt / 7)
# setting initial values and reactivity input function
ode15s.set_initial_value(u0, t=0.0).set_f_params(param)
# (DN is the number of steps between two message)
DN = 100
nsteps = DN
ii = 1
print('------------------------------')
# SOLVE THE PROBLEM
while ode15s.successful() and ode15s.t < tf:
    # this tells us how many steps have been performed
    if (ii / nsteps == 1):
        print('------------------------------')
        print(str(nsteps), ' steps performed')
        nsteps = nsteps + DN
        print('------------------------------')
    # stepwise
    ode15s.integrate(ode15s.t + dt)
    # loading
    ts.append(ode15s.t)
    sol.append(np.real(ode15s.y.copy()))
    ii = ii+1
print('------------------------------')
print('Simulation Completed')
print('------------------------------')
print('------------------------------')

######################################################################
# postprocessing #####################################################
t = np.array(ts)
sol = np.array(sol)
n = sol[:, 0]
c1 = sol[:, 1]
c2 = sol[:, 2]
c3 = sol[:, 3]
c4 = sol[:, 4]
c5 = sol[:, 5]
c6 = sol[:, 6]

plt.figure(1)
plt.plot(t, n)
plt.xlabel("Time [s]")
plt.ylabel("Population Density [1/cm]")
plt.title("Neutron concentration")
plt.grid()
plt.savefig('KINnt.png')

plt.figure(2)
plt.plot(t, c1)
plt.plot(t, c2)
plt.plot(t, c3)
plt.plot(t, c4)
plt.plot(t, c5)
plt.plot(t, c6)
plt.legend(["Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "Group 6"])
plt.xlabel("Time [s]")
plt.ylabel("Population Density [1/cm]")
plt.title("Precursor Groups concentration")
plt.grid()
plt.savefig('KINct.png')

# to show all figures
#plt.show()

# save a csv file
solFileName = "KINsol.csv"
with open(solFileName, 'w') as solFile:
  writer = csv.writer(solFile)
  rowNum = 0
  writer.writerow(['Time [s]', 'n Density [1/cm3]', 'Group 1 Precursor Density [1/cm3]  ','Group 2','Group 3','Group 4','Group 5','Group 6'])
  writer.writerows([[t[jj],n[jj],c1[jj],c2[jj],c3[jj],c4[jj],c5[jj],c6[jj]] for jj in range(len(t))])
  
