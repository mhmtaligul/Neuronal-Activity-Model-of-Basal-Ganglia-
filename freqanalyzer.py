# -*- coding: utf-8 -*-
"""
Created on Fri March  9 12:53:00 2018

@author: M.Ali GÜL
"""

from numpy import *
from brian2 import *

def freqanalyzer(data,neuronnumber,starttime,runtime):
    clippeddata = data/ms
    windwidth = 10
    stepnumber = int(((runtime-50)/windwidth) * 2)
    # Step number multiplied with 2 because next window is start to middle of before.
    freqinstep = zeros(stepnumber)
    for step in range(0,stepnumber):
        firinginstep = size(nonzero((clippeddata >= (starttime + windwidth/2*step)) & (clippeddata < (starttime + windwidth/2*(step+2)))))
        freqinstep[step] = firinginstep/neuronnumber * 1000/windwidth
    freq = mean(freqinstep)
    return freq
    
  
def spectrogram(data):
    clippeddata = data.t[find(data.t/ms>50)]/ms#First 50 ms is removed to provide stabilities
    watchedrange = 100
    spect = zeros(watchedrange) 
    for x in range(0,size(clippeddata)):
        for y in range(0,size(clippeddata)):
            rangeofspk = abs(clippeddata[x] - clippeddata[y])
            if rangeofspk != 0:
                freqofspk = int(1000 / (rangeofspk/ms))
                if freqofspk < watchedrange:
                    spect[freqofspk-1] = spect[freqofspk-1] + 1
    spect[watchedrange-1] = 0
    figure()
    plot(linspace(0,watchedrange-1,watchedrange),spect)
    return spect
    
def vectordiffer(a,b):
    return a - b 
    
def spectrogram2(data):
    data = array(data / ms, dtype = np.float32)
    periodofspks = []
    for i in range(1,size(data)):
        C = vectordiffer(data[-(size(data)-i):],data[:size(data)-i])
        periodofspks = append(periodofspks,abs(C))
    periodofspks = periodofspks[nonzero(periodofspks)]      
    frequencyofspks = 1000/periodofspks
    #frequencyofspks = detrend_mean(frequencyofspks)
    figure()
    _ = hist(frequencyofspks, 10000, histtype='step', facecolor='k', lw = 1)
    xlim(0,100)
    
def visualise_connectivity(S):
    Ns = len(S.source)
    Nt = len(S.target)
    figure(figsize=(10, 4))
    subplot(121)
    plot(zeros(Ns), arange(Ns), 'ok', ms = 5)
    plot(ones(Nt), arange(Nt), 'ok', ms = 5)
    for i, j in zip(S.i, S.j):
        if i % 2 == 1:
            plot([0, 1], [i, j], '-r')
        else:
            plot([0, 1], [i, j], '-b')
    xticks([0, 1], ['Source', 'Target'])
    ylabel('Neuron index')
    xlim(-0.1, 1.1)
    ylim(-1, max(Ns, Nt))
    
    subplot(122)
    a = [i for i in S.i]
    b = [i for i in S.j]
    plot(a, b, 'ok')
    xlim(-1, Ns)
    ylim(-1, Nt)
    xlabel('Source neuron index')
    ylabel('Target neuron index')


