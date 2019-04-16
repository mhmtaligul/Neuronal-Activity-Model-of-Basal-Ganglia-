# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 22:14:55 2018

@author: M.Ali GÜL
"""

from brian2 import *
from numpy import *
from matplotlib.pyplot import *

start_scope()
###############################################################################
#GPi model
#GPi model parameter
N_GPi = 100
C_GPi = 9
vr_GPi = -55
vt_GPi = -45
k_GPi = 1
a_GPi = 3
b_GPi = -2
c_GPi = -55
d_GPi = 150
vpeak_GPi = 20

eqsGPi = '''
dv/dt = (((k_GPi*(v-vr_GPi)*(v-vt_GPi)) - u + I)/C_GPi)/ms : 1
du/dt = a*(b*(v-vr_GPi)-u)/ms : 1
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
GGPi = NeuronGroup(N_GPi, eqsGPi, method='euler', threshold = 'v>20',reset = reset)
GGPi.vr = vr_GPi + rand(size(GGPi.a))
GGPi.vt = vt_GPi
GGPi.a = a_GPi + rand(size(GGPi.a))
GGPi.b = b_GPi + rand(size(GGPi.b)) 
GGPi.c = c_GPi + rand(size(GGPi.c))
GGPi.d = d_GPi
GGPi.v = vr_GPi
GGPi.u = GGPi.b*(GGPi.v)
###############################################################################
stm = StateMonitor(GGPi, 'v', record = True)
stm2 = StateMonitor(GGPi, 'I', record = True)
spk = SpikeMonitor(GGPi)

GGPi.I = 0
run(100*ms)

GGPi.I = 20
run(200*ms)

GGPi.I = 0
run(100*ms)

GGPi.I = 40
run(200*ms)

GGPi.I = 0
run(100*ms)

GGPi.I = -5
run(200*ms)

GGPi.I = 0
run(100*ms)
'''
GGPi.I = 100
run(1000*ms)

GGPi.I = 0
run(400*ms)
'''
figure(figsize=(18,6))
subplot(30,1,(1,24))
plot(stm.t/ms, stm.v[0], 'r', lw = 1.5)
xlim((0,1000))
ylabel('mV',fontsize = 16)
title('Globus Pallidus inside',fontsize = 32)

subplot(30,1,(26,30))
plot(stm2.t/ms, stm2.I[0], 'r', lw=3)
ylim(min(stm2.I[0]),max(stm2.I[0])+20)
xlabel('Time (ms)',fontsize = 16)
ylabel('I (uA)',fontsize = 16)

