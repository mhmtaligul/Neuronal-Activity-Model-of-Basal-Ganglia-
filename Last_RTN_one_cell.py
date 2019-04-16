# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 22:43:02 2018

@author: M.Ali GÜL
"""

from brian2 import *
from numpy import *
from matplotlib.pyplot import *

start_scope()
###############################################################################
#RTN model
N_RTN = 20
C_RTN = 5
vr_RTN = -65
vt_RTN = -45
k_RTN = 0.15
a_RTN = 0.15
b_RTN = 2
c_RTN = -55
d_RTN = 50
vpeak_RTN = 0
vb_RTN = -65

eqsRTN = '''
dv/dt = (((k_RTN*(v-vr_RTN)*(v-vt_RTN)) - u + I)/C_RTN)/ms : 1
du/dt = a*(b*(v-vr_RTN)-u)/ms : 1
vr : 1
vt : 1
a : 1
b : 1
d : 1
c : 1
I : 1
'''
reset_RTN = '''
v = c - 0.08*u
u = u + d
'''
GRTN = NeuronGroup(N_RTN, eqsRTN, method='euler', threshold = 'v>vpeak_RTN-0.08*u',reset = reset_RTN)
GRTN.vr = vr_RTN + rand(size(GRTN.vr))
GRTN.vt = vt_RTN - rand(size(GRTN.vt))
GRTN.a = a_RTN + rand(size(GRTN.a))   
GRTN.b = b_RTN + rand(size(GRTN.b))
GRTN.c = c_RTN + rand(size(GRTN.c))
GRTN.d = d_RTN - rand(size(GRTN.d))
GRTN.v = vr_RTN
GRTN.u = GRTN.b*(GRTN.v)
###############################################################################
stm = StateMonitor(GRTN, 'v', record = True)
stm2 = StateMonitor(GRTN, 'I', record = True)
spk = SpikeMonitor(GRTN)

GRTN.I = 0
run(100*ms)

GRTN.I = 50
run(200*ms)

GRTN.I = 0
run(100*ms)

GRTN.I = 100
run(200*ms)

GRTN.I = 0
run(100*ms)

GRTN.I = -100
run(200*ms)

GRTN.I = 0
run(100*ms)
'''
GRTN.I = 100
run(1000*ms)

GRTN.I = 0
run(400*ms)
'''
figure(figsize=(18,6))
subplot(30,1,(1,24))
plot(stm.t/ms, stm.v[0], 'r', lw = 1.5)
xlim((0,1000))
ylabel('mV',fontsize = 16)
title('Reticular cell',fontsize = 32)

subplot(30,1,(26,30))
plot(stm2.t/ms, stm2.I[0], 'r', lw=3)
ylim(min(stm2.I[0]),max(stm2.I[0])+20)
xlabel('Time (ms)',fontsize = 16)
ylabel('I (uA)',fontsize = 16)