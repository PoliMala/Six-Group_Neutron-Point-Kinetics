import numpy as np

# neutron kinetics data upload function definition
def dataBuilder(n0=1):
  # Problem parameter and initial condition dataBuilder
  b = 6.5e-3 * np.array([0.033, 0.219, 0.196, 0.395, 0.115, 0.042])
  l = np.array([0.0124, 0.0305, 0.111, 0.301, 1.14, 3.01])
  L = np.array([1e-5])
  c0 = b / (l * L) * n0
  u0 = np.append(np.array([n0]), c0)
  return (u0, L, l, b)

# neutron kinetics six group equations system function definition
def statedot(t, x, param):
  # Funcion: x_dot = fun(t,x,param) that returns the first derivative of
  # the system x evaluated at time t given the following param structure
  # param = [ [reactivity (callable)],
  #           [neutron mean generation time (float)],
  #           [six group dacay constant (array)],
  #           [six group resp. fractions (array)]     ]
  r = param[0]
  L = param[1]
  l = param[2]
  b = param[3]
  A = param[4]

  # computing A at time t
  A[0, :] = np.append(np.array((r(t) - b.sum()) / L), l)
  # initialize the x derivative at time t
  tmp = [0.0 for i in range(7)]  # temporary array
  x_dot = np.array(tmp)
  # computing the x derivative at time t
  for i in range(7):
    x_dot[i] = np.dot(A[i, :], np.real(x.transpose()))
  return (x_dot)


########################################################################
# Reactivity input function definition #################################

# step function
def rho_step(t, r=65e-5):
    return r  #65e-5


# single square wave function (time period DTS)
def rho_sqWave(t, DTS=10, r=650e-6):
    # reactivity function, to be defined in ppm
    if t > 0 and t < DTS:
        return r
    else:
        return 0


# periodic harmonic function
def rho_harm(t, r=650e-6):
    return r * np.sin(t)