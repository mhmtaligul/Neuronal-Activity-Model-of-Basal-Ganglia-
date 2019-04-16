# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 18:50:55 2018

@author: M.Ali GÜL
"""

from brian2 import *
from numpy import *
from matplotlib.pyplot import *

start_scope()
###############################################################################
#Posterior Cortex GABA neuron models
N_GABA = 20
C_GABA = 25
vr_GABA = -56
vt_GABA = -42
k_GABA = 0.2
a_GABA = 0.03
b_GABA = 2
c_GABA = -50
d_GABA = 2
vpeak_GABA = 40
T_GABA = 5000*ms

eqsGABA = '''
dv/dt = (((k_GABA*(v-vr_GABA)*(v-vt_GABA)) - u + I)/C_GABA)/ms : 1
du/dt = a*(b*(v-vr_GABA)-u)/ms : 1
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
GGABA = NeuronGroup(N_GABA, eqsGABA, method='euler', threshold = 'v>vpeak_GABA',reset = reset)
GGABA.vr = vr_GABA - rand(size(GGABA.vr))
GGABA.vt = vt_GABA - rand(size(GGABA.vt))
GGABA.a = a_GABA + rand(size(GGABA.a))
GGABA.b = b_GABA + rand(size(GGABA.b))
GGABA.c = c_GABA + rand(size(GGABA.c))
GGABA.d = d_GABA - rand(size(GGABA.d))
GGABA.v = -55
GGABA.I = 0
GGABA.u = GGABA.b*(GGABA.v)
###############################################################################
stm = StateMonitor(GGABA, 'v', record = True)
stm2 = StateMonitor(GGABA, 'I', record = True)
spk = SpikeMonitor(GGABA)

GGABA.I = 0
run(100*ms)

GGABA.I = 100
run(200*ms)

GGABA.I = 0
run(100*ms)

GGABA.I = 200
run(200*ms)

GGABA.I = 0
run(100*ms)

GGABA.I = -100
run(200*ms)

GGABA.I = 0
run(100*ms)
'''
GGABA.I = 100
run(1000*ms)

GGABA.I = 0
run(400*ms)
'''
figure(figsize=(18,6))
subplot(30,1,(1,24))
plot(stm.t/ms, stm.v[0], 'r', lw = 1.5)
xlim((0,1000))
ylabel('mV',fontsize = 16)
title('GABA',fontsize = 32)

subplot(30,1,(26,30))
plot(stm2.t/ms, stm2.I[0], 'r', lw=3)
ylim(min(stm2.I[0]),max(stm2.I[0])+20)
xlabel('Time (ms)',fontsize = 16)
ylabel('I (uA)',fontsize = 16)

