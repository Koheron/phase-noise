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
        
        self.n_max = 4 * self.n
        
        # Time grid (s)
        self.t = np.arange(n)/fs
        # Time steps (s)
        self.dt = 1/fs
        self.fft_freq = np.fft.fftfreq(n)/self.dt
        
        self.phase = np.zeros(self.n_max)
        self.update_phase(self.n)
        
    def update_phase(self, n_update):
        phase_steps = np.sqrt(2*self.D_phi*self.dt) * np.random.randn(n_update)
        new_phase = self.phase[-1] + np.cumsum(phase_steps)
        self.phase = np.roll(self.phase, -n_update)
        self.phase[-n_update-1:-1] = new_phase
     
    def interference(self,delay):
        return np.exp(1j * (self.phase - np.roll(self.phase, -np.int(delay / self.dt))))
        
        
        
        
fs = 125e6
n = 16384
D_phi = 1000

# Interferometer delay (s)
delay = 0.1e-6  

las = Laser(fs, n, D_phi)        

fig = plt.figure()
ax = fig.add_subplot(111)
line1, = ax.plot(np.fft.fftshift(las.fft_freq),np.linspace(-50,100,las.n))


for i in range(1000):
    las.update_phase(16384)
    signal = las.interference(delay)
    line1.set_ydata(np.fft.fftshift(10*np.log10(np.square(np.fft.fft(signal[-las.n-1:-1])))))
    fig.canvas.draw()
    time.sleep(0.01)
             

         