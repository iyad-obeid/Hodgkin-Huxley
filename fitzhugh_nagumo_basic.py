####################################
# Iyad Obeid / iyad.obeid@temple.edu
####################################

import numpy as np
import matplotlib.pyplot as plt
from neuron_models import fn_neuron ,fn_ode
from scipy.integrate import solve_ivp

# APPROACH 1 -> EXPLICITLY SOLVING ODE
def approach_1():
    n = fn_neuron()
    dt = 0.001
    t = np.arange(0,200,dt)
    v = np.zeros_like(t)
    w = np.zeros_like(t)
    for i in range(len(t)-1):
        ret_val = n.update()
        v[i+1] = ret_val[0]
        w[i+1] = ret_val[1]


    plt.figure()
    plt.plot(t,v , label='v')
    plt.plot(t,w , label='w')
    plt.title('Fitzhugh - Nagumo with ODEs solved explicitly')
    plt.legend(loc=1)
    plt.xlabel('time')


# APPROACH 2 --> USING ODE SOLVER
def approach_2():

    # define parameters
    a,b,tau,I = 0.3,1.0,12.5,0.20

    # solve ode
    soln = solve_ivp(fn_ode , [0,200] , [0,0] , args=(a,b,tau,I) )

    # plot results
    t = soln.t
    y = soln.y
    plt.figure()
    plt.plot(t,y.T)
    plt.title('Fitzhugh - Nagumo solved using SciPy ODE solver')
    plt.xlabel('time')
    plt.legend(['v','w'] , loc=1)

    plt.figure()
    plt.plot(y[0,:] , y[1,:])
    plt.xlabel('v')
    plt.ylabel('w')
    plt.title('Fitzhugh - Nagumo state space model')

if __name__ == "__main__":
    approach_1()
    approach_2()
    plt.show()