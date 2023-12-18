import sys
import os
import numpy as np
import scipy.signal
import PyQt5.QtWidgets as qtw
import PyQt5.QtCore as qtc
import PyQt5.QtGui as qtg
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui

import offset_simulation as diff

# white background, black foreground
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

class MainWindow(qtw.QMainWindow):

    def __init__(self):
        super().__init__()

        
        # Physical coordinates
        xmin = -50.0
        xmax = -xmin
        ymin = xmin
        ymax = xmax
        nx = 2000
        ny = nx

        self.chans = np.arange(16)+1     # channel number (1-16)
        self.pressure_scan_data = []
        self.efield_scan_data = []

        #self.dirnamePressure = "data/pressure_scan/20230519/"
        #self.fnamesPressure = ["05_19_2023_08_41_26_data.csv", "05_19_2023_09_31_24_data.csv",
        #                       "05_19_2023_10_17_29_data.csv", "05_19_2023_11_03_17_data.csv",
        #                       "05_19_2023_11_50_02_data.csv", "05_19_2023_12_36_24_data.csv",
        #                       "05_19_2023_13_23_33_data.csv", "05_19_2023_14_10_43_data.csv",
        #                       "05_19_2023_14_57_48_data.csv", "05_19_2023_15_45_23_data.csv",
        #                       "05_19_2023_16_33_02_data.csv", "05_19_2023_17_21_28_data.csv"
        #                       ]
        self.dirnamePressure = 'data/pressure_scan/20230519/updatedArea/'
        self.fnamesPressure = ["05_19_2023_pressureScan_200torr.csv", "05_19_2023_pressureScan_250torr.csv",
                               "05_19_2023_pressureScan_300torr.csv", "05_19_2023_pressureScan_350torr.csv",
                               "05_19_2023_pressureScan_400torr.csv", "05_19_2023_pressureScan_450torr.csv",
                               "05_19_2023_pressureScan_500torr.csv", "05_19_2023_pressureScan_550torr.csv",
                               "05_19_2023_pressureScan_600torr.csv", "05_19_2023_pressureScan_650torr.csv",
                               "05_19_2023_pressureScan_700torr.csv", "05_19_2023_pressureScan_760torr.csv" ]
        self.pressures = [200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 760]
        #
        self.dirnameEscan = "data/efield_scan/20230522/updatedArea/"
        self.fnamesEscan = ["05_22_2023_EScan_200.csv", "05_22_2023_EScan_250.csv",
                            "05_22_2023_EScan_300.csv", "05_22_2023_EScan_350.csv",
                            "05_22_2023_EScan_400.csv"]
        self.eFields = [200, 250, 300, 350, 400]

        # Configure the data plots
        self.load_pressure_scan()
        self.load_efield_scan()


        usePressureScan = False
        if usePressureScan:
            self.data_to_use = np.copy(self.pressure_scan_data)
            self.labels_to_use = self.pressures[:]
            self.nrows = 2
            self.ncols = 6
        else:
            self.data_to_use = np.copy(self.efield_scan_data)
            self.labels_to_use = self.eFields[:]
            self.nrows = 2
            self.ncols = 3
        
        self.n_data_sets = len(self.data_to_use)
        
        # Default values
        #sigx0, mux0, theta0, phi0 = 1.9, 3.3, 75, 13
        #sigdrift0 = 3.0
        sigx0, mux0, theta0, phi0 = 1.3, 3.0, 70, -25.0
        sigdrift0 = 0.3
        
        # Initialize the data
        # Image plot
        self.xx, self.yy = np.meshgrid(np.linspace(xmin, xmax, nx),
                                       np.linspace(ymin, ymax, ny))
        self.areas = diff.compute_areas()
        self.midpoints = diff.compute_midpoints()
        self.masks = diff.compute_masks2(self.xx, self.yy)

        # UI
        self.profileGraph = pg.PlotWidget()  # simulated cts. vs. radius
        self.imageDisplay = pg.PlotWidget()  # 2D image
        self.dataPlots = pg.GraphicsLayoutWidget() # 16 channels of data
        self.dataPlotAxes = []
        self.dataPlotPoints = []
        self.dataPlotLines = []
        print(f'self.n_data_sets = {self.n_data_sets}')
        #for ii in range(12): # FIXME (number of data files in pressure scan)
        for ii in range(self.n_data_sets):
            row = int(np.floor(ii/self.ncols))
            col = ii%self.ncols

            # One set of axes per dataset
            self.dataPlotAxes.append( self.dataPlots.addPlot(row=row, col=col) )
            self.dataPlotAxes[-1].setTitle(self.labels_to_use[ii])
            
            # Pressure scan or E-field scan data (points)
            self.dataPlotPoints.append( self.dataPlotAxes[-1].plot(self.chans[1:], self.data_to_use[ii][1:],
                                                                  pen=None) )
            self.dataPlotPoints[-1].setSymbol('o')

            # "fits" to data (lines)
            self.dataPlotLines.append( self.dataPlotAxes[-1].plot(self.chans[1:], self.chans[1:],
                                                                  pen=pg.mkPen(color='red', width=2) ))
            
        # Configure the 2D image plot
        self.img = pg.ImageItem()
        self.imageDisplay.addItem(self.img)
        self.img.setColorMap("viridis")
        self.imageDisplay.showAxes(True)
        self.imageDisplay.showGrid(x=True, y=True)

        # Line plot (diffusion profile)
        self.cts = np.arange(16)*10  # dummy (eventually resets per area)
        self.line = self.profileGraph.plot(self.chans, self.cts, pen=pg.mkPen(color='red', width=2))
        self.line.setSymbol('o')

        
        # Controls

        # width of fiber output illumination pattern
        self.slider_sigx_layout = qtw.QHBoxLayout()
        self.slider_sigx_label = qtw.QLabel("sigma x")
        self.sigx_factor = 10
        sigx_min, sigx_max = 0, 8  # physical units
        self.slider_sigx = qtw.QSlider(qtc.Qt.Horizontal, self)
        self.slider_sigx.setMinimum(int(sigx_min*self.sigx_factor)) # integer for slider
        self.slider_sigx.setMaximum(int(sigx_max*self.sigx_factor)) # integer for slider
        self.slider_sigx.setValue(int(sigx0*self.sigx_factor)) # integer for slider
        self.slider_sigx_layout.addWidget(self.slider_sigx_label)
        self.slider_sigx_layout.addWidget(self.slider_sigx)

        # Offset
        self.slider_mux_layout = qtw.QHBoxLayout()
        self.slider_mux_label = qtw.QLabel("mu x")
        self.mux_factor = 10
        mux_min, mux_max = 0, xmax # physical units
        self.slider_mux = qtw.QSlider(qtc.Qt.Horizontal, self)
        self.slider_mux.setMinimum(int(mux_min*self.mux_factor)) # integer for slider
        self.slider_mux.setMaximum(int(mux_max*self.mux_factor)) # integer for slider
        self.slider_mux.setValue(int(mux0*self.mux_factor)) # integer for slider
        self.slider_mux_layout.addWidget(self.slider_mux_label)
        self.slider_mux_layout.addWidget(self.slider_mux)
        
        # Theta (angle wrt cathode normal)
        self.slider_theta_layout = qtw.QHBoxLayout()
        self.slider_theta_label = qtw.QLabel("theta")
        theta_min, theta_max = 0, 85
        self.slider_theta = qtw.QSlider(qtc.Qt.Horizontal, self)
        self.slider_theta.setMinimum(theta_min)
        self.slider_theta.setMaximum(theta_max)
        self.slider_theta.setValue(int(theta0))
        self.slider_theta_layout.addWidget(self.slider_theta_label)
        self.slider_theta_layout.addWidget(self.slider_theta)

        # Phi (angle in plane of cathode)
        self.slider_phi_layout = qtw.QHBoxLayout()
        self.slider_phi_label = qtw.QLabel("phi")
        phi_min, phi_max = -90, 90
        self.slider_phi = qtw.QSlider(qtc.Qt.Horizontal, self)
        self.slider_phi.setMinimum(phi_min)
        self.slider_phi.setMaximum(phi_max)
        self.slider_phi.setValue(int(phi0))
        self.slider_phi_layout.addWidget(self.slider_phi_label)
        self.slider_phi_layout.addWidget(self.slider_phi)

        # diffusion from drift
        self.slider_sigdrift_layout = qtw.QHBoxLayout()
        self.slider_sigdrift_label = qtw.QLabel("sigma drift")
        self.sigdrift_factor = 10
        sigdrift_min, sigdrift_max = 0, 8  # physical units
        self.slider_sigdrift = qtw.QSlider(qtc.Qt.Horizontal, self)
        self.slider_sigdrift.setMinimum(int(sigdrift_min*self.sigdrift_factor)) # integer for slider
        self.slider_sigdrift.setMaximum(int(sigdrift_max*self.sigdrift_factor)) # integer for slider
        self.slider_sigdrift.setValue(int(sigdrift0*self.sigdrift_factor)) # integer for slider
        self.slider_sigdrift_layout.addWidget(self.slider_sigdrift_label)
        self.slider_sigdrift_layout.addWidget(self.slider_sigdrift)

        # checkbox to enable/disable diffusion
        self.diffusion_checkbox = qtw.QCheckBox("Apply drift")
        self.diffusion_checkbox.setChecked(True)
        
        # Layout
        self.layout = qtw.QVBoxLayout()

        self.layoutControls = qtw.QGridLayout()
        self.layoutControls.addLayout(self.slider_sigx_layout, 0, 0)
        self.layoutControls.addLayout(self.slider_mux_layout, 1, 0)
        self.layoutControls.addLayout(self.slider_theta_layout, 0, 1)
        self.layoutControls.addLayout(self.slider_phi_layout, 1, 1)
        self.layoutControls.addLayout(self.slider_sigdrift_layout, 2, 0)
        self.layoutControls.addWidget(self.diffusion_checkbox, 2, 1)
                                      
        self.layoutGraphs = qtw.QHBoxLayout()
        self.layoutGraphs.addWidget(self.imageDisplay)
        self.layoutGraphs.addWidget(self.profileGraph)

        self.layout.addLayout(self.layoutControls)
        self.layout.addLayout(self.layoutGraphs)
        self.layout.addWidget(self.dataPlots)
        
        self.setCentralWidget(qtw.QWidget())
        self.centralWidget().setLayout(self.layout)

        # Connect signals/slots
        self.slider_sigx.valueChanged.connect(self.updateGraphs)
        self.slider_mux.valueChanged.connect(self.updateGraphs)
        self.slider_theta.valueChanged.connect(self.updateGraphs)
        self.slider_phi.valueChanged.connect(self.updateGraphs)
        self.slider_sigdrift.valueChanged.connect(self.updateGraphs)
        self.diffusion_checkbox.stateChanged.connect(self.updateGraphs)

        # Display the initial plots
        self.chargeDistNoDiff = diff.make_gaussian(self.xx, self.yy, sigx0, mux0, theta0, phi0)
        self.chargeDist = np.zeros_like(self.chargeDistNoDiff)
        print("make distribution pre-diffusion")
        # apply diffusion
        self.make_diffusion_kernel(sigdrift0)
        #self.apply_drift_diffusion()
        
        self.img.setImage(self.chargeDistNoDiff.T)
        #self.img.setImage(self.kernel.T)
        self.imageDisplay.setXRange(xmin, xmax)
        self.imageDisplay.setYRange(ymin, ymax)
        self.imageDisplay.setAspectLocked(True)

        # Define the transformation (from pixel number to physical units)
        tr = qtg.QTransform()  # prepare ImageItem transformation:
        tr.scale((xmax-xmin)/nx, (ymax-ymin)/ny)       # scale horizontal and vertical axes
        tr.translate( -nx/2, -ny/2 ) # locate center at axis origin
        self.img.setTransform(tr) # assign transform

        # Default profile
        self.line.setData(x=self.chans, y=self.integrateCharge(self.chargeDist))
        self.profileGraph.setTitle(f"sigx, mux, theta, phi = {sigx0}, {mux0}, {theta0}, {phi0}")

        # overlay the anode boundaries
        #p_ellipse = pg.QtGui.QGraphicsEllipseItem(0, 0, 3, 3)  # x, y, width, height
        self.draw_anode()

        #self.add_profile_to_data_plots()
        self.updateGraphs()
        
    def add_profile_to_data_plots(self):
        xvals, yvals = self.line.getData()
        for ii in range(self.n_data_sets):
            _, ydata = self.dataPlotPoints[ii].getData()
            # scale ydata to match actual data
            self.dataPlotLines[ii].setData(xvals[1:], yvals[1:]*np.max(ydata)/np.max(yvals))
       
        
    def load_pressure_scan(self):
        #dirname = "/Users/jbattat/research/qpix/saq_analysis/data/pressure_scan/20230519/"
        fnames = [os.path.join(self.dirnamePressure, fn) for fn in self.fnamesPressure]
        #print(fnames)

        nskip = 1
        for ii in range(len(fnames)):
            junk1, rsts, junk2 = np.loadtxt(fnames[ii], delimiter=",", comments='#', skiprows=nskip, unpack=True)
            self.pressure_scan_data.append(rsts[:])

    def load_efield_scan(self):
        #dirname = "/Users/jbattat/research/qpix/saq_analysis/data/pressure_scan/20230519/"
        fnames = [os.path.join(self.dirnameEscan, fn) for fn in self.fnamesEscan]
        print(fnames)

        nskip = 1
        for ii in range(len(fnames)):
            junk1, rsts, junk2 = np.loadtxt(fnames[ii], delimiter=",", comments='#', skiprows=nskip, unpack=True)
            self.efield_scan_data.append(rsts[:])

            
    def draw_anode(self):
        for ii in range(16):
            w = 2*self.midpoints[ii]
            p_ellipse = qtw.QGraphicsEllipseItem(-w/2, -w/2, w, w)  # x, y, width, height
            p_ellipse.setPen(pg.mkPen((255, 255, 255, 100)))
            #p_ellipse.setBrush(pg.mkBrush((50, 50, 200)))
            self.imageDisplay.addItem(p_ellipse)
            
    def integrateCharge(self, cloud):
        ll = []
        for ii in range(16):
            ll.append( np.sum(cloud*self.masks[ii])/self.areas[ii])
        return ll

    def apply_drift_diffusion(self):
        #print("apply_drift_diffusion()")
        # standard convolution is ***way*** too slow
        #self.chargeDist = scipy.signal.convolve2d(self.chargeDistNoDiff, self.kernel,
        #                                          mode='full', boundary='fill', fillvalue=0)
        self.chargeDist = scipy.signal.fftconvolve(self.chargeDistNoDiff, self.kernel,
                                                   mode='same')
        
    def make_diffusion_kernel(self, sigma):
        # sigma: mm (diffusion amount from drift)
        #print('make_diffusion_kernel()')
        norm = 1.0
        self.kernel = norm*np.exp( - 0.5*(self.xx**2 + self.yy**2)/sigma**2 )
    
    def updateGraphs(self):
        sigdrift = self.slider_sigdrift.value()/self.sigdrift_factor
        sigx = self.slider_sigx.value()/self.sigx_factor
        mux = self.slider_mux.value()/self.mux_factor
        theta = self.slider_theta.value()
        phi = self.slider_phi.value()
        applyDiffusion = self.diffusion_checkbox.isChecked()
        
        #print(f'sigx, mux, theta, phi = {sigx}, {mux}, {theta}, {phi}')
        self.chargeDistNoDiff = diff.make_gaussian(self.xx, self.yy, sigx, mux, theta, phi)
        if applyDiffusion:
            self.make_diffusion_kernel(sigdrift) # update kernel
            self.apply_drift_diffusion() # convolve
        else:
            self.chargeDist = np.array([row[:] for row in self.chargeDistNoDiff])
        self.img.setImage(self.chargeDist.T)
        #self.img.setImage(self.kernel.T)
        self.line.setData(x=self.chans, y=self.integrateCharge(self.chargeDist))
        self.profileGraph.setTitle(f"sigx, mux, theta, phi, sig_drift = {sigx}, {mux}, {theta}, {phi}, {sigdrift}")

        # update lines on the pressure scan plots
        xvals, yvals = self.line.getData()
        for ii in range(self.n_data_sets):
            _, ydata = self.dataPlotPoints[ii].getData()
            # scale ydata to match actual data
            self.dataPlotLines[ii].setData(xvals[1:], yvals[1:]*np.max(ydata)/np.max(yvals))

        
if __name__ == "__main__":
    app = qtw.QApplication(sys.argv)

    win = MainWindow()
    win.show()

    sys.exit(app.exec_())
