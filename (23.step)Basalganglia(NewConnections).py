# -*- coding: utf-8 -*-
"""
Created on Mon May 29 14:47:17 2018

@author: M.Ali GÜL
"""

#For Interface
from brian2 import *
from numpy import *
from sympy import *
from matplotlib.pyplot import *
from freqanalyzer import *
from timeit import default_timer as timer
start = timer()

start_scope()
runtime = 1000 # Total simulation time which its unit will defined as ms later 
#mode is switch control variable to represent dopamine level
#mode = 0 represent the low dopamine level
#mode = 1 represent the Mediate dopamine level 
#mode = 2 represent the High dopamine level
mode = 0
DBS_switch = 1

I_dbs = 500*DBS_switch
f_dbs = 200*Hz*DBS_switch
pw_dbs = 1*ms*DBS_switch

dlevel = matrix([[0.3,0.2,50],[0.58,0.55,50],[0.70,0.9,50]])
situation = ['Low Dopamine Level','Mediate Dopamine Level','High Dopamine Level']
phi_1 = dlevel[mode,0]
phi_2 = dlevel[mode,1]
n_ext = dlevel[mode,2]

###############################################################################
#NMDA Neurons Model in Posterior Cortex
#Values of model's parameter
N_NMDA = 80
C_NMDA = 50
vr_NMDA = -60#-55
vt_NMDA = -40#-45
k_NMDA = 0.7#0.2
a_NMDA = 2#0.13
b_NMDA = -2#-1.5
c_NMDA = -60#-50
d_NMDA = 10#60
vpeak_NMDA = 35#30
fbackground_NMDA = 100 * Hz
wbackground_NMDA = 5
T_NMDA = 14000*ms

eqsNMDA = '''
dv/dt = (((k_NMDA*(v-vr_NMDA)*(v-vt_NMDA)) - u + I)/C_NMDA)/ms : 1
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
PNMDA = PoissonInput(GNMDA, 'v', 1, fbackground_NMDA, weight = wbackground_NMDA)
###############################################################################
#GABA Neurons Model in Posterior Cortex
#Values of model's parameter
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
fbackground_GABA = 100 * Hz
wbackground_GABA = 25
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
GGABA.vr = vr_GABA - rand(size(GGABA.vr))*5
GGABA.vt = vt_GABA - rand(size(GGABA.vt))*5
GGABA.a = a_GABA + rand(size(GGABA.a))*0.01
GGABA.b = b_GABA + rand(size(GGABA.b))*0.1
GGABA.c = c_GABA + rand(size(GGABA.c))*5
GGABA.d = d_GABA - rand(size(GGABA.d))
GGABA.v = -55
GGABA.I = 0
GGABA.u = GGABA.b*(GGABA.v)
###############################################################################
PGABA = PoissonInput(GGABA, 'v', 1, fbackground_GABA, wbackground_GABA)
###############################################################################
#Striatum D1 receptor's neuron
N_MSND1 = 100
C_MSND1 = 5
vr_MSND1 = -75.9
vt_MSND1 = -33.8
k_MSND1 = 1.1
a_MSND1 = 0.04
b_MSND1 = -4
c_MSND1 = -65
d_MSND1 = 0.1
vpeak_MSND1 = 40
fbackground_MSND1 = 100 * Hz
wbackground_MSND1 = 20
L_MSND1 = 0.731
gDA_MSND1 = 13.7
EDA_MSND1 = -68.4
#phi values determine the dopamine levels between 0 and 1. This level is choosen for direct pathway 
beta_MSND1 = 0.85
T_MSND1 = 10000*ms
eqsMSND1 = '''
dv/dt = (((k_MSND1*(v-vr_MSND1)*(v-vt_MSND1)) - u + I + (phi_1*gDA_MSND1*(v-EDA_MSND1)))/C_MSND1)/ms : 1
du/dt = a*(b*(v-vr_MSND1)-u)/ms : 1
vr : 1
vt : 1
a : 1
b : 1
d : 1
c : 1
I : 1
'''
reset_MSND1= '''
v = c
u = u + d
d = d*(1-L_MSND1*phi_1)
'''
GMSND1 = NeuronGroup(N_MSND1, eqsMSND1, method='euler', threshold = 'v>vpeak_MSND1',reset = reset_MSND1)
GMSND1.vr = vr_MSND1
GMSND1.vt = vt_MSND1
GMSND1.a = a_MSND1
GMSND1.b = b_MSND1
GMSND1.c = c_MSND1
GMSND1.d = d_MSND1
GMSND1.v = -65
GMSND1.u = GMSND1.b*(GMSND1.v)
###############################################################################
PMSND1 = PoissonInput(GMSND1, 'v', 1, fbackground_MSND1, weight=wbackground_MSND1)# Zamansal Poisson dağılımı ile vuru etkelyen neuron gurubu
###############################################################################
#Striatum D2 receptor's neuron
N_MSND2 = 100
C_MSND2 = 10#70
vr_MSND2 = -74#-70
vt_MSND2 = -44.1#-40
k_MSND2 = 0.9#0.25
a_MSND2 = 0.05#0.02
b_MSND2 = -15#5
c_MSND2 = -55#-65
d_MSND2 = 600#0.1
vpeak_MSND2 = 35
fbackground_MSND2 = 100 * Hz
wbackground_MSND2 = 20
alpha_MSND2 = 0.932#0.932
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
PMSND2 = PoissonInput(GMSND2, 'v', 1, fbackground_MSND2, weight=wbackground_MSND2)
###############################################################################
#Striatumda FSI interneurons
N_FSI = 30
C_FSI = 40
vr_FSI = -55#-55 
vt_FSI = -40
k_FSI = 1
a_FSI = 0.45
b_FSI = -2
c_FSI = -55
d_FSI = 2
vpeak_FSI = 25
fbackground_FSI = 100 * Hz
wbackground_FSI = 1
E_FSI = 0.525
mu_FSI = 0.65
T_FSI = 5000*ms

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
GFSI.vr = vr_FSI + rand(size(GFSI.vr))*5
GFSI.vt = vt_FSI + rand(size(GFSI.vt))*2
GFSI.a = a_FSI - rand(size(GFSI.a))*0.1 
GFSI.b = b_FSI - rand(size(GFSI.b))*0.5 
GFSI.c = c_FSI + rand(size(GFSI.c))*10 
GFSI.d = d_FSI - rand(size(GFSI.d))*0.5
GFSI.v = -65.0 + rand(size(GFSI.v))
GFSI.I = 0
GFSI.u = GFSI.b*(GFSI.v)
###############################################################################
PFSI = PoissonInput(GFSI, 'v', 1, fbackground_FSI, weight=wbackground_FSI)
###############################################################################
#GPi model
#GPi model parameters
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
fbackground_GPi = 100 * Hz
wbackground_GPi = 10
T_GPi = 20000*ms

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
GGPi.vr = vr_GPi + rand(size(GGPi.a))*10
GGPi.vt = vt_GPi
GGPi.a = a_GPi + rand(size(GGPi.a))*0.01 
GGPi.b = b_GPi + rand(size(GGPi.b)) 
GGPi.c = c_GPi + rand(size(GGPi.c))*5 
GGPi.d = d_GPi
GGPi.v = vr_GPi
GGPi.u = GGPi.b*(GGPi.v)
###############################################################################
PGPi = PoissonInput(GGPi, 'v', 1, fbackground_GPi, weight=wbackground_GPi)
###############################################################################
#GPe model
#GPe model parameters
N_GPe = 100
C_GPe = 8#5
vr_GPe = -55
vt_GPe = -45
k_GPe = 0.8#0.5
a_GPe = 0.5
b_GPe = 6
c_GPe = -50
d_GPe = 10#100 
vpeak_GPe = 25
fbackground_GPe = 100 * Hz
wbackground_GPe = 15
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
GGPe.vr = vr_GPe + rand(size(GGPe.a))*5
GGPe.vt = vt_GPe
GGPe.a = a_GPe + rand(size(GGPe.a))*0.1 
GGPe.b = b_GPe - rand(size(GGPe.b))*5 
GGPe.c = c_GPe + rand(size(GGPe.c))*10 
GGPe.v = GGPe.vr
GGPe.u = GGPe.b*(GGPe.v)
###############################################################################
PGPe = PoissonInput(GGPe, 'v', 1, fbackground_GPe, weight=wbackground_GPe)
###############################################################################
#STN modeli
#I_DBS = DBS_switch*10*(sin(2*pi*200*Hz*t)-20) : 1

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
fbackground_STN = 100 * Hz
wbackground_STN = 10
T_STN = 8000*ms

eqsSTN = '''
dv/dt = (((k_STN*(v-vr_STN)*(v-vt_STN)) - u + I + I_DBS)/C_STN)/ms : 1
du/dt = a*(b*(v-vr_STN)-u)/ms : 1
a : 1
b : 1
d : 1
c : 1
I : 1
I_DBS = DBS_switch*I_dbs*(0.5*(sign(sin(2*pi*t*f_dbs))+1))*(1-(0.5*(sign(sin(2*pi*(t+pw_dbs)*f_dbs))+1))) : 1
'''
reset = '''
v = c
u = u + d
'''
GSTN = NeuronGroup(N_STN, eqsSTN, method='euler', threshold = 'v>vpeak_STN', reset = reset)
GSTN.a = a_STN + rand(size(GSTN.a))*0.01 
GSTN.b = b_STN + rand(size(GSTN.b))*0.01 
GSTN.c = c_STN + rand(size(GSTN.c))*0.01
GSTN.d = d_STN + rand(size(GSTN.d))*0.01
GSTN.v = -85.0 + rand(size(GSTN.v))*0.01
GSTN.u = GSTN.b*(GSTN.v)
###############################################################################
PSTN = PoissonInput(GSTN, 'v', 1, fbackground_STN, weight=wbackground_STN)
###############################################################################
#Thalamus cell models
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
dv/dt = (0.04*v**2 + 5*v + 140 - u + I)/ms : 1
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
PTHL = PoissonInput(GTHL, 'v', 1, fbackground_THL, weight=wbackground_THL)
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
fbackground_RTN = 100 * Hz
wbackground_RTN = 35
T_RTN = 5000*ms

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
GRTN.vr = vr_RTN + rand(size(GRTN.vr))*10
GRTN.vt = vt_RTN - rand(size(GRTN.vt))*5
GRTN.a = a_RTN + rand(size(GRTN.a))*0.01     
GRTN.b = b_RTN + rand(size(GRTN.b))
GRTN.c = c_RTN + rand(size(GRTN.c))*10
GRTN.d = d_RTN - rand(size(GRTN.d))*10
GRTN.v = vr_RTN
GRTN.u = GRTN.b*(GRTN.v)
###############################################################################
PRTN = PoissonInput(GRTN, 'v', 1, fbackground_RTN, weight=wbackground_RTN)

syn_NMDA_NMDA = Synapses(GNMDA, GNMDA, '''w : 1
                        dge/dt = -ge / T_NMDA : 1''', 
              on_pre='''ge = ge + w*(int(ge<1))
                        v_post += ge''',method='euler')
syn_NMDA_NMDA.connect('i != j', p = 0.2)
syn_NMDA_NMDA.delay = '1*ms + 0.25*ms * randn()'
syn_NMDA_NMDA.ge = 1
syn_NMDA_NMDA.w = 0.1

syn_NMDA_GABA = Synapses(GNMDA, GGABA, '''w : 1
                        dge/dt = -ge / T_NMDA : 1''', 
              on_pre='''ge = ge + w*(int(ge<1))
                        v_post += ge''',method='euler')
syn_NMDA_GABA.connect(p = 0.2)
syn_NMDA_GABA.delay = '1*ms+rand()*ms'
syn_NMDA_GABA.ge = 1
syn_NMDA_GABA.w = 0.1

syn_NMDA_MSND1 = Synapses(GNMDA, GMSND1, '''w : 1
                        dge/dt = -ge / T_NMDA : 1''', 
              on_pre='''ge = ge + w*(int(ge<7.5))
                        v_post += ge*(1 + beta_MSND1*phi_1)''',method='euler')
syn_NMDA_MSND1.connect(p = 0.45)
syn_NMDA_MSND1.delay = '1*ms+rand()*ms'
syn_NMDA_MSND1.ge = 7.5
syn_NMDA_MSND1.w = 0.1

syn_NMDA_MSND2 = Synapses(GNMDA, GMSND2, '''w : 1
                        dge/dt = -ge / T_NMDA : 1''', 
              on_pre='''ge = ge + w*(int(ge<12.5))
                        v_post += ge*(1 - beta_MSND2*phi_2)''',method='euler')
syn_NMDA_MSND2.connect(p = 0.45)
syn_NMDA_MSND2.delay = '1*ms+rand()*ms'
syn_NMDA_MSND2.ge = 12.5
syn_NMDA_MSND2.w = 0.1

syn_NMDA_FSI = Synapses(GNMDA, GFSI, '''w : 1
                        dge/dt = -ge / T_NMDA : 1''', 
              on_pre='''ge = ge + w*(int(ge<0.1))
                        v_post += ge''',method='euler')
syn_NMDA_FSI.connect(p = 0.3)
syn_NMDA_FSI.delay = '1*ms+rand()*ms'
syn_NMDA_FSI.ge = 0.1
syn_NMDA_FSI.w = 0.1

syn_NMDA_STN = Synapses(GNMDA, GSTN, '''w : 1
                        dge/dt = -ge / T_NMDA : 1''', 
              on_pre='''ge = ge + w*(int(ge<2))
                        v_post += ge''',method='euler')
syn_NMDA_STN.connect(p = 0.2)
syn_NMDA_STN.delay = '1*ms+rand()*ms'
syn_NMDA_STN.ge = 2
syn_NMDA_STN.w = 0.1

syn_GABA_GABA = Synapses(GGABA, GGABA, '''w : 1
                        dge/dt = -ge / T_GABA : 1''', 
              on_pre='''ge = ge + w*(int(ge<0.3))
                        v_post -= ge''',method='euler')
syn_GABA_GABA.connect('i != j', p = 0.2)
syn_GABA_GABA.delay = '1*ms+rand()*ms'
syn_GABA_GABA.ge = 0.3
syn_GABA_GABA.w = 0.1

syn_GABA_NMDA = Synapses(GGABA, GNMDA, '''w : 1
                        dge/dt = -ge / T_GABA : 1''', 
              on_pre='''ge = ge + w*(int(ge<0.3))
                        v_post -= ge''',method='euler')
syn_GABA_NMDA.connect(p = 0.2)
syn_GABA_NMDA.delay = '1*ms+rand()*ms'
syn_GABA_NMDA.ge = 0.3
syn_GABA_NMDA.w = 0.1

syn_FSI_FSI = Synapses(GFSI, GFSI, '''w : 1
                        dge/dt = -ge / T_FSI : 1''', 
              on_pre='''ge = ge + w*(int(ge<0.75))
                        v_post -= ge*(0.74)''',method='euler')
syn_FSI_FSI.connect('i != j', p = 0.15)
syn_FSI_FSI.delay = '1*ms+rand()*ms'
syn_FSI_FSI.ge = 0.75
syn_FSI_FSI.w = 0.1

syn_FSI_MSND1 = Synapses(GFSI, GMSND1, '''w : 1
                        dge/dt = -ge / T_FSI : 1''', 
              on_pre='''ge = ge + w*(int(ge<0.75))
                        v_post -= ge*(1 - beta_MSND1*phi_1)''',method='euler')
syn_FSI_MSND1.connect(p = 0.15)
syn_FSI_MSND1.delay = '1*ms+rand()*ms'
syn_FSI_MSND1.ge = 0.75
syn_FSI_MSND1.w = 0.1


syn_FSI_MSND2 = Synapses(GFSI, GMSND2, '''w : 1
                        dge/dt = -ge / T_FSI : 1''', 
              on_pre='''ge = ge + w*(int(ge<0.75))
                        v_post -= ge*(1 + beta_MSND2*phi_2)''',method='euler')
syn_FSI_MSND2.connect(p = 0.15)
syn_FSI_MSND2.delay = '1*ms+rand()*ms'
syn_FSI_MSND2.ge = 0.75
syn_FSI_MSND2.w = 0.1
#?

syn_MSND1_GPi = Synapses(GMSND1, GGPi, '''w : 1
                        dge/dt = -ge / T_MSND1 : 1''', 
              on_pre='''ge = ge + w*(int(ge<25))
                        v_post -= ge''',method='euler')
syn_MSND1_GPi.connect(p = 0.3)
syn_MSND1_GPi.delay = '1*ms+rand()*ms'
syn_MSND1_GPi.ge = 25
syn_MSND1_GPi.w = 0.1

syn_MSND1_MSND1 = Synapses(GMSND1, GMSND1, '''w : 1
                        dge/dt = -ge / T_MSND1 : 1''', 
              on_pre='''ge = ge + w*(int(ge<10))
                        v_post -= ge''',method='euler')
syn_MSND1_MSND1.connect('i != j', p = 0.26)
syn_MSND1_MSND1.delay = '1*ms+rand()*ms'
syn_MSND1_MSND1.ge = 10
syn_MSND1_MSND1.w = 0.1
#?

syn_MSND2_GGPe = Synapses(GMSND2, GGPe, '''w : 1
                        dge/dt = -ge / T_MSND2 : 1''', 
              on_pre='''ge = ge + w*(int(ge<15))
                        v_post -= ge''',method='euler')
syn_MSND2_GGPe.connect(p = 0.65)
syn_MSND2_GGPe.delay = '1*ms+rand()*ms'
syn_MSND2_GGPe.ge = 15
syn_MSND2_GGPe.w = 0.1

syn_MSND2_MSND2 = Synapses(GMSND2, GMSND2, '''w : 1
                        dge/dt = -ge / T_MSND2 : 1''', 
              on_pre='''ge = ge + w*(int(ge<10))
                        v_post -= ge''',method='euler')
syn_MSND2_MSND2.connect('i != j', p = 0.15)
syn_MSND2_MSND2.delay = '1*ms+rand()*ms'
syn_MSND2_MSND2.ge = 10
syn_MSND2_MSND2.w = 0.1

syn_GPi_GPi = Synapses(GGPi, GGPi, '''w : 1
                        dge/dt = -ge / T_GPi : 1''', 
              on_pre='''ge = ge + w*(int(ge<3))
                        v_post -= ge''',method='euler')
syn_GPi_GPi.connect('i != j', p = 0.2)
syn_GPi_GPi.delay = '1*ms+rand()*ms'
syn_GPi_GPi.ge = 3
syn_GPi_GPi.w = 0.1

syn_GPi_GPe = Synapses(GGPi, GGPe, '''w : 1
                        dge/dt = -ge / T_GPi : 1''', 
              on_pre='''ge = ge + w*(int(ge<4))
                        v_post -= ge''',method='euler')
syn_GPi_GPe.connect(p = 0.4)
syn_GPi_GPe.delay = '1*ms+rand()*ms'
syn_GPi_GPe.ge = 4
syn_GPi_GPe.w = 0.1

syn_GPi_THL = Synapses(GGPi, GTHL, '''w : 1
                        dge/dt = -ge / T_GPi : 1''', 
              on_pre='''ge = ge + w*(int(ge<10))
                        v_post -= ge''',method='euler')
syn_GPi_THL.connect(p = 1)
syn_GPi_THL.delay = '1*ms+rand()*ms'
syn_GPi_THL.ge = 10
syn_GPi_THL.w = 0.1

syn_GPe_GPe = Synapses(GGPe, GGPe, '''w : 1
                        dge/dt = -ge / T_GPe : 1''', 
              on_pre='''ge = ge + w*(int(ge<2))
                        v_post -= ge''',method='euler')
syn_GPe_GPe.connect('i != j', p = 0.15)
syn_GPe_GPe.delay = '1*ms+rand()*ms'
syn_GPe_GPe.ge = 2
syn_GPe_GPe.w = 0.1

syn_GPe_STN = Synapses(GGPe, GSTN, '''w : 1
                        dge/dt = -ge / T_GPe : 1''', 
              on_pre='''ge = ge + w*(int(ge<30))
                        v_post -= ge''',method='euler')
syn_GPe_STN.connect(p = 0.2)
syn_GPe_STN.delay = '1*ms+rand()*ms'
syn_GPe_STN.ge = 30
syn_GPe_STN.w = 0.1

syn_STN_STN = Synapses(GSTN, GSTN, '''w : 1
                        dge/dt = -ge / T_STN : 1''', 
              on_pre='''ge = ge + w*(int(ge<1.5))
                        v_post += ge''',method='euler')
syn_STN_STN.connect('i !=j ', p = 0.1)
syn_STN_STN.delay = '1*ms+rand()*ms'
syn_STN_STN.ge = 1.5
syn_STN_STN.w = 0.1

syn_STN_GPe = Synapses(GSTN, GGPe, '''w : 1
                        dge/dt = -ge / T_STN : 1''', 
              on_pre='''ge = ge + w*(int(ge<2))
                        v_post += ge''',method='euler')
syn_STN_GPe.connect(p = 0.25)
syn_STN_GPe.delay = '1*ms+rand()*ms'
syn_STN_GPe.ge = 2
syn_STN_GPe.w = 0.1
#?

syn_STN_GPi = Synapses(GSTN, GGPi, '''w : 1
                        dge/dt = -ge / T_STN : 1''', 
              on_pre='''ge = ge + w*(int(ge<1))
                        v_post += ge''',method='euler')
syn_STN_GPi.connect(p = 0.3)
syn_STN_GPi.delay = '1*ms+rand()*ms'
syn_STN_GPi.ge = 1
syn_STN_GPi.w = 0.1

syn_THL_THL = Synapses(GTHL, GTHL, '''w : 1
                        dge/dt = -ge / T_THL : 1''', 
              on_pre='''ge = ge + w*(int(ge<1.5))
                        v_post += ge''',method='euler')
syn_THL_THL.connect('i != j', p = 0.2)
syn_THL_THL.delay = '1*ms+rand()*ms'
syn_THL_THL.ge = 1.5
syn_THL_THL.w = 0.1

syn_THL_RTN = Synapses(GTHL, GRTN, '''w : 1
                        dge/dt = -ge / T_THL : 1''', 
              on_pre='''ge = ge + w*(int(ge<1.5))
                        v_post += ge''',method='euler')
syn_THL_RTN.connect(p = 0.2)
syn_THL_RTN.delay = '1*ms+rand()*ms'
syn_THL_RTN.ge = 1.5
syn_THL_RTN.w = 0.1

'''
syn_THL_NMDA = Synapses(GTHL, GNMDA, 'dw/dt = -w / T_THL : 1 (event-driven)', on_pre='v_post += w')
syn_THL_NMDA.connect(p = 0.25)
syn_THL_NMDA.delay = '1*ms+rand()*ms'
syn_THL_NMDA.ge = 1
syn_THL_NMDA.w = 0
'''
syn_RTN_RTN = Synapses(GRTN, GRTN, '''w : 1
                        dge/dt = -ge / T_RTN : 1''', 
              on_pre='''ge = ge + w*(int(ge<1))
                        v_post -= ge''',method='euler')
syn_RTN_RTN.connect('i != j', p = 0.2)
syn_RTN_RTN.delay = '1*ms+rand()*ms'
syn_RTN_RTN.ge = 1
syn_RTN_RTN.w = 0.1

syn_RTN_THL = Synapses(GRTN, GTHL, '''w : 1
                        dge/dt = -ge / T_RTN : 1''', 
              on_pre='''ge = ge + w*(int(ge<1))
                        v_post -= ge''',method='euler')
syn_RTN_THL.connect(p = 0.2)
syn_RTN_THL.delay = '1*ms+rand()*ms'
syn_RTN_THL.ge = 1
syn_RTN_THL.w = 0.1

run(50*ms)

stm_NMDA = StateMonitor(GNMDA, 'v', record = True)
stm2_NMDA = StateMonitor(GNMDA, 'I', record = True)
spk_NMDA = SpikeMonitor(GNMDA)
prm_NMDA = PopulationRateMonitor(GNMDA)

stm_GABA = StateMonitor(GGABA, 'v', record = True)
stm2_GABA = StateMonitor(GGABA, 'I', record = True)
spk_GABA = SpikeMonitor(GGABA)
prm_GABA = PopulationRateMonitor(GGABA)

stm_MSND1 = StateMonitor(GMSND1, 'v', record = True)
stm2_MSND1 = StateMonitor(GMSND1, 'I', record = True)
spk_MSND1 = SpikeMonitor(GMSND1)
prm_MSND1 = PopulationRateMonitor(GMSND1)

stm_MSND2 = StateMonitor(GMSND2, 'v', record = True)
stm2_MSND2 = StateMonitor(GMSND2, 'I', record = True)
spk_MSND2 = SpikeMonitor(GMSND2)
prm_MSND2 = PopulationRateMonitor(GMSND2)

stm_FSI = StateMonitor(GFSI, 'v', record = True)
stm2_FSI = StateMonitor(GFSI, 'I', record = True)
spk_FSI = SpikeMonitor(GFSI)
prm_FSI = PopulationRateMonitor(GFSI)

stm_GPi = StateMonitor(GGPi, 'v', record = True)
stm2_GPi = StateMonitor(GGPi, 'I', record = True)
spk_GPi = SpikeMonitor(GGPi)
prm_GPi = PopulationRateMonitor(GGPi)

stm_STN = StateMonitor(GSTN, 'v', record = True)
stm2_STN = StateMonitor(GSTN, 'I', record = True)
stm_DBS_STN = StateMonitor(GSTN, 'I_DBS', record = True)
spk_STN = SpikeMonitor(GSTN)
prm_STN = PopulationRateMonitor(GSTN)

stm_GPe = StateMonitor(GGPe, 'v', record = True)
stm2_GPe = StateMonitor(GGPe, 'I', record = True)
spk_GPe = SpikeMonitor(GGPe)
prm_GPe = PopulationRateMonitor(GGPe)

stm_THL = StateMonitor(GTHL, 'v', record = True)
stm2_THL = StateMonitor(GTHL, 'I', record = True)
spk_THL = SpikeMonitor(GTHL)
prm_THL = PopulationRateMonitor(GTHL)

stm_RTN = StateMonitor(GRTN, 'v', record = True)
stm2_RTN = StateMonitor(GRTN, 'I', record = True)
spk_RTN = SpikeMonitor(GRTN)
prm_RTN = PopulationRateMonitor(GRTN)

starttime = 350 #Starttime represent the stimulation start time 
sttime = 300 #Stimulation during time

run((starttime-50)*ms) # -50*ms add because simulation already runned for 50*ms for stability
GNMDA.I += n_ext
run(sttime*ms)
GNMDA.I = 0
run((runtime-sttime-starttime)*ms)


'''
GNMDA.I += n_ext
run(runtime*ms)
'''


st_spk_NMDA = spk_NMDA.t[find(spk_NMDA.t>starttime*ms)]; st_spk_NMDA = st_spk_NMDA[find(st_spk_NMDA<=(starttime+sttime)*ms)]
st_spk_GABA = spk_GABA.t[find(spk_GABA.t>starttime*ms)]; st_spk_GABA = st_spk_GABA[find(st_spk_GABA<=(starttime+sttime)*ms)]
st_spk_MSND1 = spk_MSND1.t[find(spk_MSND1.t>starttime*ms)]; st_spk_MSND1 = st_spk_MSND1[find(st_spk_MSND1<=(starttime+sttime)*ms)]
st_spk_MSND2 = spk_MSND2.t[find(spk_MSND2.t>starttime*ms)]; st_spk_MSND2 = st_spk_MSND2[find(st_spk_MSND2<=(starttime+sttime)*ms)]
st_spk_FSI = spk_FSI.t[find(spk_FSI.t>starttime*ms)]; st_spk_FSI = st_spk_FSI[find(st_spk_FSI<=(starttime+sttime)*ms)]
st_spk_GPi = spk_GPi.t[find(spk_GPi.t>starttime*ms)]; st_spk_GPi = st_spk_GPi[find(st_spk_GPi<=(starttime+sttime)*ms)]
st_spk_STN = spk_STN.t[find(spk_STN.t>starttime*ms)]; st_spk_STN = st_spk_STN[find(st_spk_STN<=(starttime+sttime)*ms)]
st_spk_GPe = spk_GPe.t[find(spk_GPe.t>starttime*ms)]; st_spk_GPe = st_spk_GPe[find(st_spk_GPe<=(starttime+sttime)*ms)]
st_spk_THL = spk_THL.t[find(spk_THL.t>starttime*ms)]; st_spk_THL = st_spk_THL[find(st_spk_THL<=(starttime+sttime)*ms)]
st_spk_RTN = spk_RTN.t[find(spk_RTN.t>starttime*ms)]; st_spk_RTN = st_spk_RTN[find(st_spk_RTN<=(starttime+sttime)*ms)]

                         
freqNMDA = freqanalyzer(st_spk_NMDA,N_NMDA,starttime,sttime)
freqGABA = freqanalyzer(st_spk_GABA,N_GABA,starttime,sttime)
freqMSND1 = freqanalyzer(st_spk_MSND1,N_MSND1,starttime,sttime)
freqMSND2 = freqanalyzer(st_spk_MSND2,N_MSND2,starttime,sttime)
freqFSI = freqanalyzer(st_spk_FSI,N_FSI,starttime,sttime)
freqGPi = freqanalyzer(st_spk_GPi,N_GPi,starttime,sttime)
freqSTN = freqanalyzer(st_spk_STN,N_STN,starttime,sttime)
freqGPe = freqanalyzer(st_spk_GPe,N_GPe,starttime,sttime)
freqTHL = freqanalyzer(st_spk_THL,N_THL,starttime,sttime)
freqRTN = freqanalyzer(st_spk_RTN,N_RTN,starttime,sttime)

# Because of differantial equation behaviors we dont take the first 50 ms of simulation. Because the one cell model has no stability after start for a while.
clippedNMDA = spk_NMDA.t[find(spk_NMDA.t/ms>50)]/ms
clippedGABA = spk_GABA.t[find(spk_GABA.t/ms>50)]/ms
clippedMSND1 = spk_MSND1.t[find(spk_MSND1.t/ms>50)]/ms
clippedFSI = spk_FSI.t[find(spk_FSI.t/ms>50)]/ms
clippedMSND2 = spk_MSND2.t[find(spk_MSND2.t/ms>50)]/ms
clippedGPi = spk_GPi.t[find(spk_GPi.t/ms>50)]/ms
clippedSTN = spk_STN.t[find(spk_STN.t/ms>50)]/ms
clippedGPe = spk_GPe.t[find(spk_GPe.t/ms>50)]/ms
clippedTHL = spk_THL.t[find(spk_THL.t/ms>50)]/ms
clippedRTN = spk_RTN.t[find(spk_RTN.t/ms>50)]/ms

histstep = 100


def firingprint(spk,clippeddata,N_group,freqgroup,tit=""):
    plot(spk.t/ms, spk.i, '.k',label = int(freqgroup)*Hz, ms = 8)
    xlim(0,runtime)
    ylabel('Nerve Cell',fontsize = 36,)
    title(tit,fontsize = 48)
    legend(loc = 1,fontsize = 36)
    [hist1,hist2,hist3] = hist(spk.t/ms, histstep, histtype='step', facecolor='b', lw = 4, label = 'Spiking rate')
    legend(loc = 1,fontsize = 36)
    savefig('Figures/'+tit+'.png')
    return hist1
    
def firingprint2(spk,clippeddata,N_group,freqgroup,tit=""):
    plot(spk.t/ms, spk.i, '.k',label = int(freqgroup)*Hz, ms = 2)
    xlim(0,runtime)
    ylabel('Nerve Cell',fontsize = 12)
    title(tit,fontsize = 12)
    legend(loc = 1,fontsize = 12)
    [hist1,hist2,hist3] = hist(spk.t/ms, histstep, histtype='step', facecolor='b', lw = 1.5, label = 'Spiking rate')
    legend(loc = 1,fontsize = 12)    
    return hist1

def allprint():
    rc('xtick', labelsize=8); rc('ytick', labelsize=8)
    
    figure(figsize=(20,10))
    subplot(4,3,1); hist_NMDA = firingbastır2(spk_NMDA,clippedNMDA,N_NMDA,freqNMDA,tit="NMDA Cells")
    subplot(4,3,2); hist_GABA = firingbastır2(spk_GABA,clippedGABA,N_GABA,freqGABA,tit="GABA Cells")
    subplot(4,3,4); hist_MSND1 = firingbastır2(spk_MSND1,clippedMSND1,N_MSND1,freqMSND1,tit="D1 type optimal")
    subplot(4,3,5); hist_FSI = firingbastır2(spk_FSI,clippedFSI,N_FSI,freqFSI,tit="Fast Spiking Rates")
    subplot(4,3,6); hist_MSND2 = firingbastır2(spk_MSND2,clippedMSND2,N_MSND2,freqMSND2,tit="D2 type optimal")
    subplot(4,3,7); hist_GPi = firingbastır2(spk_GPi,clippedGPi,N_GPi,freqGPi,tit="Globus Pallidus inside")
    subplot(4,3,8); hist_STN = firingbastır2(spk_STN,clippedSTN,N_STN,freqSTN,tit="Subtalamic Cell")
    subplot(4,3,9); hist_GPe = firingbastır2(spk_GPe,clippedGPe,N_GPe,freqGPe,tit="Globus Pallidus outside")
    subplot(4,3,10); hist_THL = firingbastır2(spk_THL,clippedTHL,N_THL,freqTHL,tit="Thalamus")
    subplot(4,3,11); hist_RTN = firingbastır2(spk_RTN,clippedRTN,N_RTN,freqRTN,tit="Reticular Cell")
    subplot(4,3,12); plot(stm_DBS_STN.t/ms,stm_DBS_STN.I_DBS[0],'k', lw = 2, label = "DBS")
    text(130,(I_dbs + 100)*0.2,'Freq = '+str(f_dbs)+"\n"+'Spiking time = '+str(pw_dbs),fontsize=10,bbox={'facecolor':'white', 'alpha':1, 'pad':10})
    title("DBS spikes",fontsize = 12)
    xlim(100,150)
    ylim(-10,I_dbs + 100)
    legend(loc = 1,fontsize = 12)
    
    subplot(4,3,3); plot(stm2_NMDA.t/ms,stm2_NMDA.I[0],'r', lw = 2, label = "Applied notice")
    title("Cortex Spikes notice",fontsize = 12)
    ylim(0,max(stm2_NMDA.I[0])+10)
    legend(loc = 1,fontsize = 12)
       
    suptitle(situation[mode],fontsize = 24)
    #savefig('Figures/NewConnections'+str(counter)+'.trial.png')
    '''
    figure(figsize=(20,10))
    
    subplot(4,3,1); psd(hist_NMDA,65536,histstep)
    title('Frequency analysis via LPF'); show()
    
    
    subplot(4,3,2); psd(hist_GABA,65536,histstep)
    title('Frequency analysis via LPF'); show()
    
    subplot(4,3,4); psd(hist_MSND1,65536,histstep)
    title('Frequency analysis via LPF'); show()
    
    subplot(4,3,5); psd(hist_FSI,65536,histstep)
    title('Frequency analysis via LPF'); show()
    
    subplot(4,3,6); psd(hist_MSND2,65536,histstep)
    title('Frequency analysis via LPF'); show()
    
    subplot(4,3,7); psd(hist_GPi,65536,histstep)
    title('Frequency analysis via LPF'); show()
    
    subplot(4,3,8); psd(hist_STN,65536,histstep)
    title('Frequency analysis via LPF'); show()
    
    subplot(4,3,9); psd(hist_GPe,65536,histstep)
    title('Frequency analysis via LPF'); show()
    
    subplot(4,3,10); psd(hist_THL,65536,histstep)
    title('Frequency analysis via LPF'); show()
    
    subplot(4,3,11); psd(hist_RTN,65536,histstep)
    title('Frequency analysis via LPF'); show()
    suptitle(situation[mode]) 
    '''
def differprint():
    rc('xtick', labelsize=20); rc('ytick', labelsize=20)#Grafik eksenlerinde yazan değerlerin fon büyüklüğünü belirler.

    figure(figsize=(18,6)); hist_NMDA = firingbastır(spk_NMDA,clippedNMDA,N_NMDA,freqNMDA,tit="NMDA Cells")
    figure(figsize=(18,6)); hist_GABA = firingbastır(spk_GABA,clippedGABA,N_GABA,freqGABA,tit="GABA Cells")
    figure(figsize=(18,6)); hist_MSND1 = firingbastır(spk_MSND1,clippedMSND1,N_MSND1,freqMSND1,tit="D1 type optimal")
    figure(figsize=(18,6)); hist_FSI = firingbastır(spk_FSI,clippedFSI,N_FSI,freqFSI,tit="Fast Spiking Rates")
    figure(figsize=(18,6)); hist_MSND2 = firingbastır(spk_MSND2,clippedMSND2,N_MSND2,freqMSND2,tit="D2 type optimal")
    figure(figsize=(18,6)); hist_GPi = firingbastır(spk_GPi,clippedGPi,N_GPi,freqGPi,tit="Globus Pallidus inside")
    figure(figsize=(18,6)); hist_STN = firingbastır(spk_STN,clippedSTN,N_STN,freqSTN,tit="Subtalamic Celll")
    figure(figsize=(18,6)); hist_GPe = firingbastır(spk_GPe,clippedGPe,N_GPe,freqGPe,tit="Globus Pallidus outside")
    figure(figsize=(18,6)); hist_THL = firingbastır(spk_THL,clippedTHL,N_THL,freqTHL,tit="Thalamus")
    figure(figsize=(18,6)); hist_RTN = firingbastır(spk_RTN,clippedRTN,N_RTN,freqRTN,tit="Reticular Cell")

    figure(figsize=(18,3.2))
    plot(stm_DBS_STN.t/ms,stm_DBS_STN.I_DBS[0],'k', lw = 2, label = "DBS")
    text(143,(I_dbs + 100)*0.2,'Freq = '+str(f_dbs)+"\n"+'Spiking time'= '+str(pw_dbs),fontsize=16,bbox={'facecolor':'white', 'alpha':1, 'pad':10})
    title("DBS notice",fontsize = 24)
    xlim(100,150)
    ylim(-10,I_dbs + 100)
    legend(loc = 1,fontsize = 24)
    savefig('Figures/DBS.png')

    figure(figsize=(18,2.7))
    plot(stm2_NMDA.t/ms,stm2_NMDA.I[0],'r', lw = 4, label = "Applied notice")
    ylim(0,max(stm2_NMDA.I[0])+10)
    legend(loc = 1,fontsize = 24)
    savefig('Figures/notice.png')
    

allprint()

#visualise_connectivity(syn_NMDA_MSND1) #sdfdfasdf

basalganglia_time = timer() - start

print("Simulation Time = %f" % basalganglia_time)