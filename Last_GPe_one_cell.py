# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 22:23:35 2018

@author: M.Ali GÜL
"""

from brian2 import *
from numpy import *
from matplotlib.pyplot import *

start_scope()
###############################################################################
#GPe model
N_GPe = 100
C_GPe = 8
vr_GPe = -55
vt_GPe = -45
k_GPe = 0.8
a_GPe = 0.5
b_GPe = 6
c_GPe = -50
d_GPe = 10
vpeak_GPe = 25
T_GPe = 1000*ms

eqsGPe = '''
dv/dt = (((k_GPe*(v-vr_GPe)*(v-vt_GPe)) - u + I)/C_GPe)/ms : 1
du/dt = a*(b*(v-vr_GPe)-u)/ms : 1
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
GGPe = NeuronGroup(N_GPe, eqsGPe, method='euler', threshold = 'v>vpeak_GPe',reset = reset)
GGPe.vr = vr_GPe + rand(size(GGPe.a))
GGPe.vt = vt_GPe
GGPe.a = a_GPe + rand(size(GGPe.a))
GGPe.b = b_GPe - rand(size(GGPe.b))
GGPe.c = c_GPe + rand(size(GGPe.c))
GGPe.d = d_GPe - rand(size(GGPe.d))
GGPe.v = GGPe.vr
GGPe.u = GGPe.b*(GGPe.v)
###############################################################################
stm = StateMonitor(GGPe, 'v', record = True)
stm2 = StateMonitor(GGPe, 'I', record = True)
spk = SpikeMonitor(GGPe)

GGPe.I = 0
run(100*ms)

GGPe.I = 100
run(200*ms)

GGPe.I = 0
run(100*ms)

GGPe.I = 200
run(200*ms)

GGPe.I = 0
run(100*ms)

GGPe.I = -100
run(200*ms)

GGPe.I = 0
run(100*ms)
'''
GGPe.I = 100
run(1000*ms)

GGPe.I = 0
run(400*ms)
'''
figure(figsize=(18,6))
subplot(30,1,(1,24))
plot(stm.t/ms, stm.v[0], 'r', lw = 1.5)
xlim((0,1000))
ylabel('mV',fontsize = 16)
title('Globus Pallidus outside',fontsize = 32)

subplot(30,1,(26,30))
plot(stm2.t/ms, stm2.I[0], 'r', lw=3)
ylim(min(stm2.I[0]),max(stm2.I[0])+20)
xlabel('Time (ms)',fontsize = 16)
ylabel('I (uA)',fontsize = 16)

