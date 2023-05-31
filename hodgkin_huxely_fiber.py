####################################
# Iyad Obeid / iyad.obeid@temple.edu
####################################

import numpy as np
import matplotlib.pyplot as plt
from neuron_models import hh_neuron


# DEFINE TIME #
dt   = 0.01
tmax = 30.0 
t = np.arange(0.0, tmax, dt)

# DEFINE NUM OF NEURONS # 
nCells = 60
my_neurons = [ hh_neuron(dt) for i in range(nCells) ]

# DEFINE CONDUCTION BETWEEN NEURONS
gConnect = 5

# DEFINE VECTORS #
v     = np.zeros((nCells,len(t)))
IStim = np.zeros(len(t))

# DEFINE STIM STRENGTH AND DURATION #
tSTIM_START = 5
tSTIM_DUR = 5
STIM_STRENGTH = 20
IStim[ (t>tSTIM_START) & (t<=tSTIM_START+tSTIM_DUR) ] = STIM_STRENGTH

# DEFINE MISC FUNCTIONS
def isFirstCell(j): return j == 0
def isLastCell(j):  return j == (nCells-1)

# RUN SIMULATION
for ith_time in range(len(t)-1): 

    if ith_time%100 == 0: print(t[ith_time])

    for ith_neuron,n in enumerate(my_neurons):

        # DETERMINE STIMULUS CURRENT 
        Istimulus = 0 
        if isFirstCell( ith_neuron ):
            Istimulus = IStim[ ith_time ]
        
        # DETERMINE CURRENT FROM CELL TO THE LEFT #
        I_Left = 0
        if not isFirstCell(ith_neuron):
            I_Left = gConnect * (v[ith_neuron-1, ith_time] - v[ith_neuron, ith_time])
        
        # DETERMINE CURRENT FROM CELL TO THE RIGHT #
        I_Right = 0
        if not isLastCell(ith_neuron):  
            I_Right = gConnect * (v[ith_neuron+1, ith_time] - v[ith_neuron,ith_time])

        # UPDATE NEURON     
        v[ith_neuron,ith_time+1] = n.update(Istimulus+I_Left+I_Right)



# PLOT IT #
f,a = plt.subplots(2,1,sharex=True)
a[0].plot(t,v[0:nCells:5,:].T) 
a[0].set_ylabel('membrane voltage (mv)')
a[0].set_title('membrane voltage, every fifth neuron')
a[1].plot(t,IStim)
a[1].set_xlabel('time (ms)')
a[1].set_ylabel('stimulus current (mA)')
plt.show()
