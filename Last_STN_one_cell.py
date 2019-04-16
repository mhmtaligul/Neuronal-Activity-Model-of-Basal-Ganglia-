# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 22:26:11 2018

@author: M.Ali GÜL
"""

from brian2 import *
from numpy import *
from matplotlib.pyplot import *

start_scope()
###############################################################################
#STN model
N_STN = 100
C_STN = 20
vr_STN = -60
vt_STN = -45
k_STN = 1
a_STN = 0.005
b_STN = 0.265
c_STN = -65
d_STN = 2
vpeak_STN = 30
vb_STN = -65

eqsSTN = '''
dv/dt = (((k_STN*(v-vr_STN)*(v-vt_STN)) - u + I)/C_STN)/ms : 1
du/dt = a*(b*(v-vr_STN)-u)/ms : 1
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
GSTN = NeuronGroup(N_STN, eqsSTN, method='euler', threshold = 'v>vpeak_STN', reset = reset)
GSTN.a = a_STN + rand(size(GSTN.a))
GSTN.b = b_STN + rand(size(GSTN.b))
GSTN.c = c_STN + rand(size(GSTN.c))
GSTN.d = d_STN + rand(size(GSTN.d))
GSTN.v = -85.0 + rand(size(GSTN.v))
GSTN.u = GSTN.b*(GSTN.v)
###############################################################################
stm = StateMonitor(GSTN, 'v', record = True)
stm2 = StateMonitor(GSTN, 'I', record = True)
spk = SpikeMonitor(GSTN)

GSTN.I = 0
run(100*ms)

GSTN.I = 100
run(200*ms)

GSTN.I = 0
run(100*ms)

GSTN.I = 200
run(200*ms)

GSTN.I = 0
run(100*ms)

GSTN.I = -100
run(200*ms)

GSTN.I = 0
run(100*ms)
'''
GSTN.I = 100
run(1000*ms)

GSTN.I = 0
run(400*ms)
'''
figure(figsize=(18,6))
subplot(30,1,(1,24))
plot(stm.t/ms, stm.v[0], 'r', lw = 1.5)
xlim((0,1000))
ylabel('mV',fontsize = 16)
title('Subthalamic Cell',fontsize = 32)

subplot(30,1,(26,30))
plot(stm2.t/ms, stm2.I[0], 'r', lw=3)
ylim(min(stm2.I[0]),max(stm2.I[0])+20)
xlabel('Time (ms)',fontsize = 16)
ylabel('I (uA)',fontsize = 16)

