#!/usr/bin/env python3

import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtWidgets


class Laser:
    def __init__(self, fs, n, D_phi):
        self.fs = float(fs)
        self.dt = 1 / self.fs
        self.n = n
        self.D_phi = D_phi  # Phase diffusion coefficient (rad^2/s).

        # Keep twice the FFT length so the plotted samples have a delay history.
        self.n_ = 2 * n
        self.t_ = np.arange(self.n_) * self.dt
        self.f_fft_ = np.fft.fftfreq(self.n_, self.dt)
        self.f_fft = np.fft.fftfreq(n, self.dt)
        self.phase = np.full(self.n_, 2 * np.pi * np.random.rand())
        self.update_phase(self.n_)
        self.interferometer_phase = np.pi / 4

    def update_phase(self, n_update):
        # dphi = sqrt(2 D_phi dt) N(0, 1); optical FWHM = D_phi / pi.
        steps = np.sqrt(2 * self.D_phi * self.dt) * np.random.randn(n_update)
        last_phase = self.phase[-1]
        self.phase = np.roll(self.phase, -n_update)
        self.phase[-n_update:] = last_phase + np.cumsum(steps)

    def shift_phase(self, delay, correct_slope=False):
        # FFT shifts assume a periodic record. Remove its drift before shifting.
        phase = self.phase
        if correct_slope:
            endpoint = phase[-1] + (phase[1] - phase[0] + phase[-1] - phase[-2]) / 2
            slope = (endpoint - phase[0]) / (self.n_ * self.dt)
            phase = phase - slope * self.t_
        # exp(-i 2 pi f tau) gives phi(t - tau).
        spectrum = np.fft.fft(phase) * np.exp(-2j * np.pi * self.f_fft_ * delay)
        shifted = np.fft.ifft(spectrum).real
        if correct_slope:
            shifted += slope * (self.t_ - delay)
        return shifted

    def interference_signal(self, delay):
        # Equal-power arms: I(t) = [1 + cos(bias + phi(t) - phi(t - tau))] / 2.
        phase_difference = self.phase - self.shift_phase(delay, correct_slope=True)
        return 0.5 * (1 + np.cos(self.interferometer_phase + phase_difference))


class Window(QtWidgets.QWidget):
    def __init__(self, laser):
        super().__init__()
        self.laser = laser
        self.setWindowTitle("Koheron Simulation of laser phase noise")
        self.resize(800, 600)

        self.dphi_step = np.pi * 1e3  # 1 kHz optical FWHM per slider step.
        self.dphi_label = QtWidgets.QLabel()
        self.dphi_slider = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
        self.dphi_slider.setRange(0, int(10e6 / self.dphi_step))
        self.dphi_slider.setValue(round(laser.D_phi / self.dphi_step))

        self.delay_step = 0.1e-9  # 0.1 ns per slider step, up to 100 ns.
        self.delay_label = QtWidgets.QLabel()
        self.delay_slider = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
        self.delay_slider.setRange(0, round(100e-9 / self.delay_step))
        self.delay_slider.setValue(self.delay_slider.maximum())

        self.plot = pg.PlotWidget()
        self.plot.setLabel("bottom", "Frequency (MHz)")
        self.plot.setLabel("left", "Intensity PSD (dB/Hz)")
        self.plot.setMouseEnabled(x=False, y=True)
        self.curve = self.plot.plot(pen=(0, 4))

        controls = QtWidgets.QFormLayout()
        controls.addRow(self.dphi_label, self.dphi_slider)
        controls.addRow(self.delay_label, self.delay_slider)
        layout = QtWidgets.QVBoxLayout(self)
        layout.addLayout(controls)
        layout.addWidget(self.plot)

        self.dphi_slider.valueChanged.connect(self.change_dphi)
        self.delay_slider.valueChanged.connect(self.change_delay)
        self.change_dphi()
        self.change_delay()
        self.update_plot()
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.update_plot)
        self.timer.start(30)

    def change_dphi(self):
        self.laser.D_phi = self.dphi_slider.value() * self.dphi_step
        self.dphi_label.setText(f"Linewidth FWHM (kHz): {self.laser.D_phi / (np.pi * 1e3):.2f}")

    def change_delay(self):
        self.delay = self.delay_slider.value() * self.delay_step
        self.delay_label.setText(f"Delay (ns): {1e9 * self.delay:.2f}")

    def update_plot(self):
        self.laser.update_phase(16)
        signal = self.laser.interference_signal(self.delay)[-self.laser.n:]
        window = np.hanning(self.laser.n)
        spectrum = np.fft.fft((signal - signal.mean()) * window)
        # Two-sided PSD: its integral over frequency is the windowed mean square.
        psd = np.abs(spectrum)**2 / (self.laser.fs * np.sum(window**2))
        psd_db = 10 * np.log10(np.maximum(psd, np.finfo(float).tiny))
        self.curve.setData(1e-6 * np.fft.fftshift(self.laser.f_fft), np.fft.fftshift(psd_db))


def main():
    fs = 125e6  # Sampling frequency (Hz).
    n = 1024  # Samples per spectrum.
    linewidth = 100e3  # Optical FWHM (Hz).
    app = pg.mkQApp()
    app.setQuitOnLastWindowClosed(True)
    window = Window(Laser(fs, n, D_phi=np.pi * linewidth))
    window.show()
    return app.exec()


if __name__ == '__main__':
    import sys
    sys.exit(main())
