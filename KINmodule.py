import numpy as np

# neutron kinetics data upload function definition
def dataBuilder(n0=1.0,b0=6.5e-3,fj=[0.033,0.219,0.196,0.395,0.115,0.042],\
                lj=[0.0124,0.0305,0.111,0.301,1.14,3.01],L0=1e-5):
  # Problem parameter and initial condition dataBuilder
  b = b0 * np.array(fj)
  l = np.array(lj)
  L = np.array([L0])
  c0 = b / (l * L) * n0
  u0 = np.append(np.array([n0]), c0)
  return (u0, L, l, b)

# neutron kinetics six group equations system function definition
def KINdot(t, x, param):
  # Funcion: x_dot = fun(t,x,param) that returns the first derivative of
  # the system x evaluated at time t given the following param structure
  # param = [ [reactivity (callable)],
  #           [neutron mean generation time (float)],
  #           [six group dacay constant (array)],
  #           [six group resp. fractions (array)]     ]
  r = param[0]
  L = param[1]
  b = param[2]
  A = param[3]
  q = param[4]
  p = param[-1]
  # computing A at time t
  A[0][0] = (r[p[0]] - b.sum())/L
  # initialize the x derivative at time t
  xdot = np.array([0.0 for ii in range(7)])
  # computing the x derivative at time t
  for i in range(7):
    xdot[i] = np.dot(A[i, :], np.real(x.transpose()))
  xdot[1] = xdot[1]+q
  return (xdot)


########################################################################
# Reactivity input function definition #################################

# step function
def rho_step(t, r=65e-5, t0=0):
    if t>t0:
        return r  #65e-5
    else:
        return 0

# ramp function
def rho_ramp(t, rt=65e-6, t0=0):
    if t>t0:
        return rt*(t-t0)
    else:
        return 0

# single square wave function (time period DTS)
def rho_sqWave(t, r=650e-6, t0=0, DTS=10):
    # reactivity function, to be defined in ppm
    if t > t0 and t < DTS:
        return r
    else:
        return 0


# periodic harmonic function
def rho_harm(t, r=650e-6, t0=0, T=1):
    return r * np.sin((t-t0)*2*np.pi/T)
