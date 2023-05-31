####################################
# Iyad Obeid / iyad.obeid@temple.edu
####################################

import numpy as np

class hh_neuron:

    def __init__(self,dt):
        self.Ek  , self.gk  = -12 , 36
        self.Ena , self.gna = 115 , 120
        self.El  , self.gl  = 10.6 , 0.3
        self.C   = 1

        self.v   = 0
        self.m   = self.alpha_m() / (self.alpha_m() + self.beta_m() )
        self.h   = self.alpha_h() / (self.alpha_h() + self.beta_h() )
        self.n   = self.alpha_n() / (self.alpha_n() + self.beta_n() )

        self.dt = dt

    def alpha_m(self):
        return (2.5-0.1*self.v)/(np.exp(2.5-0.1*self.v)-1)
    def beta_m(self):
        return 4*np.exp(-self.v/18)
    def alpha_h(self):
        return 0.07*np.exp(-self.v/20)
    def beta_h(self):
        return 1/(np.exp(3-0.1*self.v)+1)
    def alpha_n(self):
        return (0.1-0.01*self.v)/(np.exp(1-0.1*self.v)-1)
    def beta_n(self): 
        return 0.125*np.exp(-self.v/80)

    def update(self , I_stim):
        I_na = self.gna * (self.m**3) * self.h * (self.v-self.Ena)
        I_k  = self.gk  * (self.n**4)          * (self.v-self.Ek)
        I_l  = self.gl                         * (self.v-self.El)
        I_mem = I_na + I_k + I_l

        dv = (I_stim - I_mem) / self.C
        dm = self.alpha_m() * (1-self.m) - self.beta_m()*self.m
        dh = self.alpha_h() * (1-self.h) - self.beta_h()*self.h
        dn = self.alpha_n() * (1-self.n) - self.beta_n()*self.n

        self.v += dv * self.dt
        self.m += dm * self.dt
        self.h += dh * self.dt
        self.n += dn * self.dt

        return self.v

class fn_neuron:

    def __init__(self , a=0.3 , b=1.0 , tau=12.5 , I=0.20 , dt=0.001):
        self.v , self.w = 0,0
        self.a , self.b = a,b
        self.tau = tau
        self.I   = I
        self.dt  = dt

    def update(self):
        v,w = self.v , self.w
        a,b = self.a , self.b
        tau = self.tau
        I = self.I
        dt = self.dt

        dv = v - (v**3)/3 - w + I
        dw = (v + a -b*w)/tau

        self.v += dv*dt
        self.w += dw*dt

        return self.v , self.w

def fn_ode(t,y,a,b,tau,I):
    v,w = y   
    dv  = v - (v**3)/3 - w + I
    dw  = (v + a - b*w)/tau

    return np.array( [dv,dw] )
