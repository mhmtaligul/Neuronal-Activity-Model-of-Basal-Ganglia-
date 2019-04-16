# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 22:37:45 2018

@author: M.Ali GÜL
"""

from brian2 import *
from numpy import *
from matplotlib.pyplot import *

start_scope()
###############################################################################
#Thalamic Cell Model
N_THL = 80
C_THL = 13
vr_THL = -60
vt_THL = -50
k_THL = 0.11
a_THL = 0.05
b_THL = 2
c_THL = -65
d_THL = 2
vpeak_THL = 35

eqsTHL = '''
dv/dt = (((k_THL*(v-vr_THL)*(v-vt_THL)) - u + I)/C_THL)/ms : 1
du/dt = a*(b*(v-vr_THL)-u)/ms : 1
vr : 1
vt : 1
a : 1
b : 1
d : 1
c : 1
I : 1
'''
reset_THL = '''
v = c - 0.1*u
u = u + d
'''
GTHL = NeuronGroup(N_THL, eqsTHL, method='euler', threshold = 'v>vpeak_THL+0.1*u',reset = reset_THL)
GTHL.vr = vr_THL
GTHL.vt = vt_THL
GTHL.a = a_THL
GTHL.b = b_THL
GTHL.c = c_THL
GTHL.d = d_THL
GTHL.v = vr_THL
GTHL.u = GTHL.b*(GTHL.v)
###############################################################################
stm = StateMonitor(GTHL, 'v', record = True)
stm2 = StateMonitor(GTHL, 'I', record = True)
spk = SpikeMonitor(GTHL)

GTHL.I = 0
run(100*ms)

GTHL.I = 100
run(200*ms)

GTHL.I = 0
run(100*ms)

GTHL.I = 200
run(200*ms)

GTHL.I = 0
run(100*ms)

GTHL.I = -100
run(200*ms)

GTHL.I = 0
run(100*ms)
'''
GTHL.I = 100
run(1000*ms)

GTHL.I = 0
run(400*ms)
'''
figure(figsize=(18,6))
subplot(30,1,(1,24))
plot(stm.t/ms, stm.v[0], 'r', lw = 1.5)
xlim((0,1000))
ylabel('mV',fontsize = 16)
title('Thalamus',fontsize = 32)

subplot(30,1,(26,30))
plot(stm2.t/ms, stm2.I[0], 'r', lw=3)
ylim(min(stm2.I[0]),max(stm2.I[0])+20)
xlabel('Time (ms)',fontsize = 16)
ylabel('I (uA)',fontsize = 16)