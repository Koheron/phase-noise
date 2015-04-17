#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import time

# https://www.institutoptique.fr/content/download/3234/22015/file/Optique%20Statistique%20cours%20ecrit.pdf

class Laser:
        
    def __init__(self, fs, n, D_phi):    
        # Sampling frequency (Hz)
        self.fs = fs
        # Number of points in time    
        self.n = n

        self.D_phi = D_phi   
        
        self.n_max = 2 * self.n
        
        # Time grid (s)
        self.t = np.arange(n)/fs
        # Time steps (s)
        self.dt = 1/fs
        self.fft_freq = np.fft.fftfreq(n)/self.dt
        
        self.phase = np.zeros(self.n_max) + 2* np.pi * np.random.rand()
        self.update_phase(self.n)

        
    def update_phase(self, n_update):
        phase_steps = np.sqrt(2*self.D_phi*self.dt) * np.random.randn(n_update)
        self.phase = np.roll(self.phase, -n_update)
        self.phase[-n_update-1:-1] = np.cumsum(phase_steps) + self.phase[-n_update-2]
     
    def interference_signal(self,delay):
        return np.real(np.exp(1j * (self.phase - np.roll(self.phase, np.int(delay / self.dt)))))
        
        
        
        
fs = 125e6
n = 16384
D_phi = 100000

# Interferometer delay (s)
delay = 0.1e-6  

las = Laser(fs, n, D_phi)        

fig = plt.figure()
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)
line1, = ax1.plot(np.fft.fftshift(las.fft_freq),np.linspace(-50,100,las.n))
line2, = ax2.plot(las.t,np.linspace(0.9,1.1,las.n))
line3, = ax3.plot(las.t,np.linspace(-10.1,10.1,las.n))



for i in range(100):
    las.update_phase(128)
    signal = las.interference_signal(delay)
    line1.set_ydata(np.fft.fftshift(10*np.log10(np.square(np.fft.fft(signal[-las.n-1:-1])))))
    line2.set_ydata(signal[-las.n-1:-1])   
    line3.set_ydata(las.phase[-las.n-1:-1])
    fig.canvas.draw()
    time.sleep(0.01)
             

         