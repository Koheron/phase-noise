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
        
        self.interferometer_phase = np.pi/4

        
    def update_phase(self, n_update):
        phase_steps = np.sqrt(2*self.D_phi*self.dt) * np.random.randn(n_update)
        self.phase = np.roll(self.phase, -n_update)
        self.phase[-n_update-1:-1] = np.cumsum(phase_steps) + self.phase[-n_update-2]
     
    def interference_signal(self,delay):
        signal = np.real(np.exp(1j * self.interferometer_phase) * \
                 np.exp(1j * (self.phase - np.roll(self.phase, np.int(delay / self.dt))))) 
        return signal
        
fs = 250e6
n = 1024
linewidth = 100e3; 
D_phi = 2*np.pi*linewidth

# Interferometer delay (s)
delay = 0.1e-6  

las = Laser(fs, n, D_phi)        

# Plot settings
fig = plt.figure(figsize=(16,12))
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)
line1, = ax1.plot(np.fft.fftshift(las.fft_freq),np.linspace(-50,50,las.n))
line11, = ax1.plot(np.fft.fftshift(las.fft_freq),np.linspace(-50,50,las.n))
line2, = ax2.plot(las.t,np.linspace(-1.1,1.1,las.n))
line3, = ax3.plot(las.t,np.linspace(-10.1,10.1,las.n))
line31, = ax3.plot(las.t,np.linspace(-10.1,10.1,las.n))

psd_avg = np.zeros(las.n)

for i in range(1000):
    las.update_phase(128)
    signal = las.interference_signal(delay)
    psd = np.abs(np.square(np.fft.fft(signal[-las.n-1:-1])));
    psd_avg = (i * psd_avg + psd)/(i+1) 
    line1.set_ydata(np.fft.fftshift(10*np.log10(psd)))
    line11.set_ydata(np.fft.fftshift(10*np.log10(psd_avg)))
    line2.set_ydata(signal[-las.n-1:-1])   
    line3.set_ydata(las.phase[-las.n-1:-1])
    tmp = np.roll(las.phase, np.int(delay / las.dt))
    line31.set_ydata(tmp[-las.n-1:-1])
    fig.canvas.draw()
    time.sleep(0.01)
             

         