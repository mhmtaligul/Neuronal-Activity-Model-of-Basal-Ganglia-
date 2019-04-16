# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 21:43:13 2018

@author: M.Ali GÜL
"""

from brian2 import *
from numpy import *
from matplotlib.pyplot import *

start_scope()
###############################################################################
#Striatum FSI interneurons
N_FSI = 30
C_FSI = 40
vr_FSI = -55
vt_FSI = -40
k_FSI = 1
a_FSI = 0.45
b_FSI = -2
c_FSI = -55
d_FSI = 2
vpeak_FSI = 25
E_FSI = 0.525
mu_FSI = 0.65

reset = '''
v = c
u = u + d
'''
eqsFSI = '''
dv/dt = (((k_FSI*(v-(vr_FSI*0.72))*(v-vt_FSI)) - u + I)/C_FSI)/ms : 1
du/dt = a*(b*(v-vr_FSI)-u)/ms : 1
vr : 1
vt : 1
a : 1
b : 1
d : 1
c : 1
I : 1
'''

GFSI = NeuronGroup(N_FSI, eqsFSI, method='euler', threshold = 'v>vpeak_FSI',reset = reset)
GFSI.vr = vr_FSI + rand(size(GFSI.vr))
GFSI.vt = vt_FSI + rand(size(GFSI.vt))
GFSI.a = a_FSI - rand(size(GFSI.a))
GFSI.b = b_FSI - rand(size(GFSI.b))
GFSI.c = c_FSI + rand(size(GFSI.c))
GFSI.d = d_FSI - rand(size(GFSI.d))
GFSI.v = -65.0 + rand(size(GFSI.v))
GFSI.I = 20
GFSI.u = GFSI.b*(GFSI.v)
###############################################################################
stm = StateMonitor(GFSI, 'v', record = True)
stm2 = StateMonitor(GFSI, 'I', record = True)
spk = SpikeMonitor(GFSI)

GFSI.I = 1
run(100*ms)

GFSI.I = 100
run(200*ms)

GFSI.I = 0
run(100*ms)

GFSI.I = 200
run(200*ms)

GFSI.I = 0
run(100*ms)

GFSI.I = -50
run(200*ms)

GFSI.I = 0
run(100*ms)
'''
GFSI.I = 100
run(1000*ms)

GFSI.I = 0
run(400*ms)
'''
figure(figsize=(18,6))
subplot(30,1,(1,24))
plot(stm.t/ms, stm.v[0], 'r', lw = 1.5)
xlim((0,1000))
ylabel('mV',fontsize = 16)
title('Fast Spiking Interneurons',fontsize = 32)

subplot(30,1,(26,30))
plot(stm2.t/ms, stm2.I[0], 'r', lw=3)
ylim(min(stm2.I[0]),max(stm2.I[0])+20)
xlabel('Time (ms)',fontsize = 16)
ylabel('I (uA)',fontsize = 16)

