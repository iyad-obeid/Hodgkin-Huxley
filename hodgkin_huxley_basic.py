####################################
# Iyad Obeid / iyad.obeid@temple.edu
####################################

import numpy as np
import matplotlib.pyplot as plt
from neuron_models import hh_neuron


# DEFINE TIME #
dt   = 0.01
tmax = 30.0 
t    = np.arange(0.0, tmax, dt)

# DEFINE VECTORS #
v     = np.zeros(len(t))
Istim = np.zeros(len(t))

# DEFINE STIM STRENGTH AND DURATION #
tSTIM_START   =  5
tSTIM_DUR     =  5
STIM_STRENGTH = 20
Istim[ (t>tSTIM_START) & (t<=tSTIM_START+tSTIM_DUR) ] = STIM_STRENGTH

# CREATE A NEURON
n = hh_neuron(dt)

# SIMULATE THE NEURON
for i in range(len(t)-1): 
   
    if i%100 == 0: 
        print(t[i])

    v[i+1] = n.update(Istim[i])


# PLOT IT #
f,a = plt.subplots(2,1,sharex=True)
a[0].plot(t,v)
a[0].set_ylabel('membrane voltage (mv)')
a[1].plot(t,Istim)
a[1].set_xlabel('time (s)')
a[1].set_ylabel('stimulus current (mA)')

plt.show()
