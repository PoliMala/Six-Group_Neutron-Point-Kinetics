#!/usr/bin/python3
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
import csv
import sys
import getopt as gop
from KINmodule import *
import os
#  set the 'MPLCONFIGDIR' to '/tmp' to get a writable environment for matplotlib
os.environ['MPLCONFIGDIR'] = '/tmp'

def main(argv):
    ############################## DEFAULT OPTIONS ############################
    # data
    [u0, L, l, b] = dataBuilder()
    # time discretizzation
    t0 = 0.0
    tf = 100.0
    dt = 0.005
    # source strenght
    q = 0
    # reactivity function
    rho_fun = rho_step
    # plots
    PLOT = True
    # output Files
    OUT = False
    try: ############ CHECK THE CORRECT PROGRAM USAGE #########################
        opts, args = gop.getopt(argv,'d:t:p:o:q:r:h',['data','time','plot',\
        'out','source','reac'])
    except gop.GetoptError:
            print('usage: main.py -i <inputfile>')
            sys.exit(1)
    ########################### ARGUMENT PARSING ##############################
    for opt, arg in opts:
        if opt == '-h':
    ######################## TIME DISCRETIZZATION #############################
            print('usage: python -opt <value> main.py\n\
the options are ...\n\n')
            sys.exit(0)
    ############################## DATA LOADING ###############################
        elif opt in ('-d','--data'):
            if arg == "TRIGA":
                [u0, L, l, b] = dataBuilder(1.0,7.3e-3)
            else:
                [u0, L, l, b] = dataBuilder()
    ############################## PLOT OPTION ################################
        elif opt in ('-p','--plot'):
            PLOT = (arg=='1')
    ########################### SOURCE TERM OPTION ############################
        elif opt in ('-o','--out'):
            q = float(arg)
    ########################### OUTPUT CSV OPTION #############################
        elif opt in ('-q','--source'):
            OUT = (arg=='1')
    ##################### REACTIVITY FUNCTION SELECTION #######################
        elif opt in ('-r','--reac'):
            ropt = arg.split()
            if ropt[0] == 'step':
                ti = float(ropt[1])
                r0 = float(ropt[2])
                def rho_fun(t):
                    return rho_step(t,ti,r0)
            if ropt[0] == 'ramp':
                ti = float(ropt[1])
                rt = float(ropt[2])
                def rho_fun(t):
                    return rho_ramp(t,ti,rt)
            elif ropt[0] == 'harm':
                ti = float(ropt[1])
                T  = float(ropt[2])
                r0 = float(ropt[3])
                def rho_fun(t):
                    return rho_harm(t,ti,T,r0)
            elif ropt[0] == 'sqW':
                ti = float(ropt[1])
                DT = float(ropt[2])
                r0 = float(ropt[3])
                def rho_fun(t):
                    return rho_sqWave(t,ti,DT,r0)
    ######################## TIME DISCRETIZZATION #############################
        elif opt in ('-t','--time'):
            topt = arg.split()
            t0 = float(topt[0])
            tf = float(topt[1])
            dt = float(topt[2])
    ########################### INITIALIZZATION ###############################
    N   = round((tf-t0)/dt)+1
    t = np.array(range(N))*dt+t0
    # time index pointer
    p = np.array([0])
    # reactivity function init
    r = np.array([0.0 for ii in range(N)])
    # initializing matrix A: KINdot system matrix (input dependent)
    tmp = [[0.0 for i in range(7)] for j in range(7)]  # temporary array
    A = np.array(tmp)
    # constant part of matrix A
    A[0, :] = np.append(np.array(0.0), l)
    for i in range(6):
        A[(i + 1), 0] = b[i] / L
        A[(i + 1), (i + 1)] = -l[i]
    # set the KINdot parameters (see KINmodule)
    param = [r, L, b, A, q, p]
    ######################## PROBLEM RESOLUTION ###############################
    # SOLUTION initializzation
    sol = [u0]
    ## SOLUTOR parameters (same solutor as MATLAB ode15s )
    # System matrix
    solver = ode(KINdot)
    # max step limitation necessary
    solver.set_integrator('lsoda',max_step = dt/7)
    # setting initial values and reactivity input function
    solver.set_initial_value(u0, t=t0).set_f_params(param)
    # (DN is the number of steps between two message)
    DN = int(N/5)
    nsteps = DN
    print('===================================================================')
    # SOLVE THE PROBLEM
    for ii in range(0,N-1):
        # this tells us how many steps have been performed
        if (ii / nsteps == 1):
            print('===================================================================')
            print(str(nsteps), ' steps performed out of ',str(N-1))
            nsteps = nsteps + DN
            print('===================================================================')
        if solver.successful() and solver.t < tf:
            p[0] = ii+1 # begin from time-step 1 (time-step 0 -> initial cond.)
            # loading the reactivity
            r[p[0]] = rho_fun(t[p[0]])
            # stepwise
            solver.integrate(t[p[0]])
            # loading the solution
            sol.append(np.real(solver.y.copy()))
    if PLOT:
        # postprocessing #####################################################
        sol = np.array(sol)
        # neutron population and precursors concentration extraction
        n = sol[:, 0]
        c1 = sol[:, 1]
        c2 = sol[:, 2]
        c3 = sol[:, 3]
        c4 = sol[:, 4]
        c5 = sol[:, 5]
        c6 = sol[:, 6]
        # Neutron population plot
        plt.figure(1)
        plt.plot(t, n)
        plt.xlabel("Time [s]")
        plt.ylabel("Population Density [1/cm]")
        plt.title("Neutron concentration")
        plt.grid()
        plt.savefig('../KINnt.png')
        # Delayed group plot
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
        plt.savefig('../KINct.png')
        # to show all figures
        #plt.show()
    # save a csv file
    solFileName = "../KINsol.csv"
    with open(solFileName, 'w') as solFile:
      writer = csv.writer(solFile)
      rowNum = 0
      writer.writerow(['Time [s]', 'n Density [1/cm3]', 'Group 1 Precursor Density [1/cm3]  ','Group 2','Group 3','Group 4','Group 5','Group 6'])
      writer.writerows([[t[jj],n[jj],c1[jj],c2[jj],c3[jj],c4[jj],c5[jj],c6[jj]] for jj in range(len(t))])

# Execute if not imported
if __name__ == '__main__':
    main(sys.argv[1:])
    print('===================================================================')
    print('Simulation Completed')
    print('===================================================================')
    print('===================================================================')