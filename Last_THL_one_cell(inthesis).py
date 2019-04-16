# -*- coding: utf-8 -*-
"""
Created on Tue May  8 02:08:00 2018

@author: M.Ali GÜL
"""

from brian2 import *
from numpy import *
from matplotlib.pyplot import *

start_scope()
###############################################################################
#Thalamic Cell Model
N_THL = 80
a_THL = 0.05
b_THL = 0.26
c_THL = -60
d_THL = 0
vpeak_THL = 35
fbackground_THL = 100 * Hz
wbackground_THL = 10
T_THL = 6000*ms

eqsTHL = '''
dv/dt = (0.04*v**2 + 5*v + 140 - u + I)/ms + 0.5 * (xi*ms**-0.5) : 1
du/dt = a*((b*v) - u)/ms : 1
a : 1
b : 1
c : 1
d : 1
I : 1
'''
reset_THL = '''
v = c
u = u + d
'''
GTHL = NeuronGroup(N_THL, eqsTHL, method='euler', threshold = 'v>vpeak_THL+0.1*u',reset = reset_THL)
GTHL.a = a_THL
GTHL.b = b_THL 
GTHL.c = c_THL 
GTHL.d = d_THL 
GTHL.v = -65.0 
GTHL.I = 0
GTHL.u = GTHL.b*(GTHL.v)
###############################################################################

stm = StateMonitor(GTHL, 'v', record = True)
stm2 = StateMonitor(GTHL, 'I', record = True)
spk = SpikeMonitor(GTHL)

GTHL.I = 0
run(100*ms)

GTHL.I = 1
run(200*ms)

GTHL.I = 0
run(100*ms)

GTHL.I = 2
run(200*ms)

GTHL.I = 0
run(100*ms)

GTHL.I = -10
run(50*ms)

GTHL.I = 0
run(100*ms)

GTHL.I = -10
run(50*ms)

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
