# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 21:31:29 2018

@author: M.Ali GÜL
"""

from brian2 import *
from numpy import *
from matplotlib.pyplot import *
#0.2 0.55 0.9
phi_2 = 0.2

start_scope()
###############################################################################
#Striatum D2 receptors neuron
N_MSND2 = 100
C_MSND2 = 10
vr_MSND2 = -74
vt_MSND2 = -44.1
k_MSND2 = 0.9
a_MSND2 = 0.05
b_MSND2 = -15
c_MSND2 = -55
d_MSND2 = 600
vpeak_MSND2 = 35
alpha_MSND2 = 0.932
#phi values determine the dopamine levels between 0 and 1. This level is choosen for direct pathway 
beta_MSND2 = 0.156
E_DA2 = -68.0
gDA2 = 31.1
T_MSND2 = 10000*ms
eqsMSND2 = '''
dv/dt = (((k_MSND2*(v-vr_MSND2)*(v-vt_MSND2)) - u + I + (phi_2*gDA2*(v-E_DA2)))/C_MSND2)/ms : 1
du/dt = a*(b*(v-vr_MSND2)-u)/ms : 1
vr : 1
vt : 1
a : 1
b : 1
d : 1
c : 1
I : 1
'''
reset = '''
v = c
u = u + d
'''

GMSND2 = NeuronGroup(N_MSND2, eqsMSND2, method='euler', threshold = 'v>vpeak_MSND2',reset = reset)
GMSND2.vr = vr_MSND2
GMSND2.vt = vt_MSND2
GMSND2.a = a_MSND2
GMSND2.b = b_MSND2
GMSND2.c = c_MSND2
GMSND2.d = d_MSND2
GMSND2.v = -65.0
GMSND2.u = GMSND2.b*(GMSND2.v)
###############################################################################
stm = StateMonitor(GMSND2, 'v', record = True)
stm2 = StateMonitor(GMSND2, 'I', record = True)
spk = SpikeMonitor(GMSND2)

GMSND2.I = 0
run(100*ms)

GMSND2.I = 150
run(200*ms)

GMSND2.I = 0
run(100*ms)

GMSND2.I = 335
run(200*ms)

GMSND2.I = 0
run(100*ms)

GMSND2.I = -100
run(200*ms)

GMSND2.I = 0
run(100*ms)
'''
GMSND2.I = 100
run(1000*ms)

GMSND2.I = 0
run(400*ms)
'''
figure(figsize=(18,6))
subplot(30,1,(1,24))
plot(stm.t/ms, stm.v[0], 'r', lw = 1.5)
text(900,0,'phi = '+str(phi_2),fontsize=16,bbox={'facecolor':'white', 'alpha':1, 'pad':10})
xlim((0,1000))
ylabel('mV',fontsize = 16)
title('D2 type optimal cell',fontsize = 32)

subplot(30,1,(26,30))
plot(stm2.t/ms, stm2.I[0], 'r', lw=3)
ylim(min(stm2.I[0]),max(stm2.I[0])+20)
xlabel('Time (ms)',fontsize = 16)
ylabel('I (uA)',fontsize = 16)

