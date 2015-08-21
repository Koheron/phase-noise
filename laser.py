#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg

# https://www.institutoptique.fr/content/download/3234/22015/file/Optique%20Statistique%20cours%20ecrit.pdf

class Laser:
        
    def __init__(self, fs, n, D_phi):    
        # Sampling frequency (Hz)
        self.fs = fs
        # Number of points in time    
        self.n = n

        self.D_phi = D_phi   
        
        self.n_ = 2 * self.n
        
        # Time grid (s)
        self.t  = np.arange(self.n)/fs
        self.t_ = np.arange(self.n_)/fs
        # Time steps (s)
        self.dt = 1/fs
        # Frequency grid (Hz)
        self.f_fft = np.fft.fftfreq(self.n)/self.dt        
        self.f_fft_ = np.fft.fftfreq(self.n_)/self.dt
        
        self.phase = np.zeros(self.n_) + 2* np.pi * np.random.rand()
        self.update_phase(self.n)
        
        self.interferometer_phase = np.pi/4

    def update_phase(self, n_update):
        phase_steps = np.sqrt(2*self.D_phi*self.dt) * np.random.randn(n_update)
        self.phase = np.roll(self.phase, -n_update)
        self.phase[-n_update-1:-1] = np.cumsum(phase_steps) + self.phase[-n_update-2]
     
    def interference_signal(self, delay):
        signal = 0.5*(1 + np.real(np.exp(1j * self.interferometer_phase) * \
                 np.exp(1j * (self.phase - self.shift_phase(delay, correct_slope=True)))))
        return signal
        
    def shift_phase(self, delay, correct_slope=False):
        phase = self.phase
        wedge = np.exp(2*1j*np.pi*delay*self.f_fft_)
        if correct_slope:
            data_inter = phase[-1] + (phase[1]-phase[0]+phase[-1]-phase[-2])/2
            slope = (data_inter - phase[0]) / (self.n_ * self.dt)
            phase = phase - slope * self.t_
        fft_phase = np.fft.fft(phase)
        phase_shifted = np.real(np.fft.ifft(fft_phase * wedge))
        if correct_slope:
            phase_shifted += slope * self.t_
        return phase_shifted
        

class Window(QtGui.QMainWindow):
    def __init__(self, app, laser):
        super(Window, self).__init__()
        self.app = app        
        self.laser = laser
        
        self.max_dphi = 10e6
        self.dphi_step = 10e3
        
        self.delay = 0.1e-6
        self.max_delay = 0.1e-6
        self.delay_step = 0.1e-9
                
        self.setWindowTitle("Koheron Simulation of laser phase noise") # Title
        self.setWindowIcon(QtGui.QIcon('icon_koheron.png'))
        self.resize(800, 600) # Size

        # Layout
        self.centralWid = QtGui.QWidget()
        self.setCentralWidget(self.centralWid)
        
        self.lay = QtGui.QVBoxLayout()
        self.button_sublayout = QtGui.QHBoxLayout()   
        self.hlay1 = QtGui.QHBoxLayout()   

        self.value_layout = QtGui.QVBoxLayout()         
        self.slider_layout = QtGui.QVBoxLayout()                
        
        self.centralWid.setLayout(self.lay)
                
        # Widgets : Buttons, Sliders, PlotWidgets
      
        # D_phi
        self.dphi_label = QtGui.QLabel()
        self.dphi_label.setText('Linewitdth (kHz): '+"{:.2f}".format(self.laser.D_phi/(2*np.pi)))
        self.dphi_slider = QtGui.QSlider()
        self.dphi_slider.setMinimum(0)
        self.dphi_slider.setMaximum(self.max_dphi/self.dphi_step)
        self.dphi_slider.setOrientation(QtCore.Qt.Horizontal)
        
        # D_phi
        self.delay_label = QtGui.QLabel()
        self.delay_label.setText('Delay (ns): '+"{:.2f}".format(self.delay))
        self.delay_slider = QtGui.QSlider()
        self.delay_slider.setMinimum(0)
        self.delay_slider.setMaximum(self.max_delay/self.delay_step)
        self.delay_slider.setOrientation(QtCore.Qt.Horizontal)        
        
        
        # Plot Widget   
        self.plotWid = pg.PlotWidget(name="data")
        self.dataItem = pg.PlotDataItem(1e-6 * np.fft.fftshift(self.laser.f_fft),0*self.laser.t, pen=(0,4))
        self.plotWid.addItem(self.dataItem)
        self.plotItem = self.plotWid.getPlotItem()
        self.plotItem.setMouseEnabled(x=False, y = True)
        #specItem.setYRange(-8192, 8192)
        # Axis
        self.plotAxis = self.plotItem.getAxis("bottom")
        self.plotAxis.setLabel("Frequency (MHz)")
        
        # Add Widgets to layout
        
        self.value_layout.addWidget(self.dphi_label,0)        
        self.slider_layout.addWidget(self.dphi_slider,0)
        
        self.value_layout.addWidget(self.delay_label,0)        
        self.slider_layout.addWidget(self.delay_slider,0)        
        
        
        self.hlay1.addLayout(self.value_layout)     
        self.hlay1.addLayout(self.slider_layout) 
        
        self.lay.addLayout(self.hlay1)
        self.lay.addWidget(self.plotWid)
        
        self.dphi_slider.valueChanged.connect(self.change_dphi)
        self.delay_slider.valueChanged.connect(self.change_delay)
        
        self.show()
        
        # Define events
        
    def update(self):        
        self.laser.update_phase(16)
        signal = self.laser.interference_signal(self.delay)
        psd = np.abs(np.square(np.fft.fft(signal[-self.laser.n-1:-1])));
        self.dataItem.setData(1e-6 * np.fft.fftshift(self.laser.f_fft),np.fft.fftshift(10*np.log10(psd)))
    
  
    def change_dphi(self):
        self.laser.D_phi = self.dphi_slider.value()*self.dphi_step    
        self.dphi_label.setText('Linewidth (kHz) : '+"{:.2f}".format(1e-3 * self.laser.D_phi / (2*np.pi)))
        
    def change_delay(self):
        self.delay = self.delay_slider.value()*self.delay_step 
        self.delay_label.setText('Delay (ns) : '+"{:.2f}".format(1e9 * self.delay))    
        
        
def main():        
    fs = 125e6
    n = 1024
    linewidth = 100e3; # Laser linewidth (Hz)
    D_phi = 2*np.pi*linewidth    
   
    las = Laser(fs, n, D_phi)        
    
    app = QtGui.QApplication.instance()
    if app == None:
        app = QtGui.QApplication([])
    app.quitOnLastWindowClosed()    
    
    win = Window(app, las)    
    
    while True:         
        win.update()
        QtGui.QApplication.processEvents()

if __name__ == '__main__':
    import sys    
    main()
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()     

         