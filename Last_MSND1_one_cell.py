# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 18:52:42 2018

@author: M.Ali GÜL
"""
from brian2 import *
from numpy import *
from matplotlib.pyplot import *
#0.3 0.58 0.7
phi_1 = 0.7

start_scope()
###############################################################################
#Striatum D1 receptors neuron
N_MSND1 = 100
C_MSND1 = 5
vr_MSND1 = -75.9
vt_MSND1 = -33.8
k_MSND1 = 1.1
a_MSND1 = 0.04
b_MSND1 = -4
c_MSND1 = -65
d_MSND1 = 0.1
vpeak_MSND1 = 40
L_MSND1 = 0.731
gDA_MSND1 = 13.7
EDA_MSND1 = -68.4
#phi values determine the dopamine levels between 0 and 1. This level is choosen for direct pathway 
beta_MSND1 = 0.85
T_MSND1 = 10000*ms
eqsMSND1 = '''
dv/dt = (((k_MSND1*(v-vr_MSND1)*(v-vt_MSND1)) - u + I + (phi_1*gDA_MSND1*(v-EDA_MSND1)))/C_MSND1)/ms : 1
du/dt = a*(b*(v-vr_MSND1)-u)/ms : 1
vr : 1
vt : 1
a : 1
b : 1
d : 1
c : 1
I : 1
'''
reset_MSND1= '''
v = c
u = u + d
d = d*(1-L_MSND1*phi_1)
'''
GMSND1 = NeuronGroup(N_MSND1, eqsMSND1, method='euler', threshold = 'v>vpeak_MSND1',reset = reset_MSND1,dt = 0.1*ms)
GMSND1.vr = vr_MSND1
GMSND1.vt = vt_MSND1
GMSND1.a = a_MSND1
GMSND1.b = b_MSND1
GMSND1.c = c_MSND1
GMSND1.d = d_MSND1
GMSND1.v = -65
GMSND1.u = GMSND1.b*(GMSND1.v)
###############################################################################
stm = StateMonitor(GMSND1, 'v', record = True)
stm2 = StateMonitor(GMSND1, 'I', record = True)
spk = SpikeMonitor(GMSND1)

GMSND1.I = 0
run(100*ms)

GMSND1.I = 150
run(200*ms)

GMSND1.I = 0
run(100*ms)

GMSND1.I = 335
run(200*ms)

GMSND1.I = 0
run(100*ms)

GMSND1.I = -100
run(200*ms)

GMSND1.I = 0
run(100*ms)
'''
GMSND1.I = 100
run(1000*ms)

GMSND1.I = 0
run(400*ms)
'''
figure(figsize=(18,6))
subplot(30,1,(1,24))
plot(stm.t/ms, stm.v[0], 'r', lw = 1.5)
text(900,20,'phi = '+str(phi_1),fontsize=16,bbox={'facecolor':'white', 'alpha':1, 'pad':10})
xlim((0,1000))
ylabel('mV',fontsize = 16)
title('D1 type optimal cell',fontsize = 32)

subplot(30,1,(26,30))
plot(stm2.t/ms, stm2.I[0], 'r', lw=3)
ylim(min(stm2.I[0]),max(stm2.I[0])+20)
xlabel('Time (ms)',fontsize = 16)
ylabel('I (uA)',fontsize = 16)

