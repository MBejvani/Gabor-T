# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 19:11:39 2017

@author: Bj
"""
from numpy import *
import scipy,pylab
#save('stft_sig',s,dt)
s = load('stft_sig.npy');dt = load('dt.npy');
s = s.T;   s = s[0];   n = len(s)
t = dt*linspace(0,n-1,n)
f = linspace(0,1/(2*dt),n//2)
fn = 1/(2*dt)
df = f[4]-f[3]
sd_max = max(t)/2/pi/df
sd_min = dt**2/pi
def AB_STF(s,dt,sd):
    n = len(s)
    t = dt*linspace(0,n-1,n)
    TF = []
    for i in range(n):
        if shape(sd)==():
            sdd = sd
        else:
            sdd = sd[i]
        w = 1/(sqrt(2*pi)*sdd) * exp(-1/2*((t-t[i])/sdd)**2)
        TF.append(scipy.fft(s*w))
    ATF = absolute(TF)
    ATF = ATF[:,n//2:]
    return  ATF

sdl = arange(sd_min,sd_max,sd_max/30)
CF = []
for sd in sdl:
    ATF = AB_STF(s,dt,sd)
    ATF = log(ATF+0.001)
    CFi = []
    m=20
    for i in range(n):
        a = i-m
        b = i+m
        if i<m:
            a = 0
        if i>n-m:
            b = n
        CFi.append(sum(ATF[a:b,:]))
    CF.append(CFi)
CF = array(CF)
sd_opt = []
sd_in_fig = []
for i in range(n):
    CFy = CF[:,i]
    CFy = list(CFy)
    sd_in = CFy.index(min(CFy))
    sd_in_fig.append(sd_in)
    sd_opt.append(sdl[sd_in])

TT = AB_STF(s,dt,sd_opt) #0.002*ones([n,1])
#TT = array(TT).transpose() # axis r wrong - esp==>for instantaneous
pylab.figure()
pylab.imshow(TT, origin='lower', aspect='auto',
             interpolation='nearest')
pylab.figure()
pylab.subplot(121)
pylab.imshow(CF, origin='lower', aspect='auto',interpolation='nearest')
pylab.subplot(122)
pylab.plot(sd_in_fig)
pylab.show()
# sd=1/2pi/f'_ins ==> max_sd = t_max/2pi/df & min_sd = dt/2pi/fn
#sd_opt,sd_max,sd_min,mean([sd_max,sd_min])