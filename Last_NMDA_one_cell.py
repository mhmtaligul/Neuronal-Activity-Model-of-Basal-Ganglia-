# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 12:53:00 2018

@author: M.Ali GÜL
"""

from brian2 import *
from numpy import *
from matplotlib.pyplot import *

start_scope()
###############################################################################
#Posterior Cortex NMDA neuron models
N_NMDA = 80
C_NMDA = 50
vr_NMDA = -60
vt_NMDA = -40
k_NMDA = 0.7
a_NMDA = 2
b_NMDA = -2
c_NMDA = -60
d_NMDA = 10
vpeak_NMDA = 35

eqsNMDA = '''
dv/dt = (((k_NMDA*(v-vr)*(v-vt)) - u + I)/C_NMDA)/ms : 1
du/dt = a*(b*(v-vr_NMDA)-u)/ms : 1
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

GNMDA = NeuronGroup(N_NMDA, eqsNMDA, method='euler', threshold = 'v>vpeak_NMDA',reset = reset)
GNMDA.vr = vr_NMDA + rand(size(GNMDA.vr))
GNMDA.vt = vt_NMDA + rand(size(GNMDA.vt))
GNMDA.a = a_NMDA + rand(size(GNMDA.a))
GNMDA.b = b_NMDA + rand(size(GNMDA.b))
GNMDA.c = c_NMDA + rand(size(GNMDA.c))
GNMDA.d = d_NMDA - rand(size(GNMDA.d))
GNMDA.v = -65
GNMDA.u = GNMDA.b*(GNMDA.v)
###############################################################################
stm = StateMonitor(GNMDA, 'v', record = True)
stm2 = StateMonitor(GNMDA, 'I', record = True)
spk = SpikeMonitor(GNMDA)

GNMDA.I = 0
run(100*ms)

GNMDA.I = 100
run(200*ms)

GNMDA.I = 0
run(100*ms)

GNMDA.I = 200
run(200*ms)

GNMDA.I = 0
run(100*ms)

GNMDA.I = -100
run(200*ms)

GNMDA.I = 0
run(100*ms)
'''
GNMDA.I = 100
run(1000*ms)

GNMDA.I = 0
run(400*ms)
'''
figure(figsize=(18,6))
subplot(30,1,(1,24))
plot(stm.t/ms, stm.v[0], 'r', lw = 1.5)
xlim((0,1000))
ylabel('mV',fontsize = 16)
title('NMDA',fontsize = 32)

subplot(30,1,(26,30))
plot(stm2.t/ms, stm2.I[0], 'r', lw=3)
ylim(min(stm2.I[0]),max(stm2.I[0])+20)
xlabel('Time (ms)',fontsize = 16)
ylabel('I (uA)',fontsize = 16)

