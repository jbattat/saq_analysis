import sys
import os
import numpy as np
import PyQt5.QtWidgets as qtw
import PyQt5.QtCore as qtc
import PyQt5.QtGui as qtg
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui

import offset_simulation as diff

class MainWindow(qtw.QMainWindow):

    def __init__(self):
        super().__init__()

        # white background, black foreground
        pg.setConfigOption('background', 'w')
        pg.setConfigOption('foreground', 'k')
        
        # Physical coordinates
        xmin = -50.0
        xmax = -xmin
        ymin = xmin
        ymax = xmax
        nx = 2000
        ny = nx

        self.chans = np.arange(16)+1     # channel number (1-16)
        self.pressure_scan_data = []

        self.dirname = "data/pressure_scan/20230519/"
        self.fnames = ["05_19_2023_08_41_26_data.csv", "05_19_2023_11_50_02_data.csv",
                       "05_19_2023_14_57_48_data.csv", "05_19_2023_09_31_24_data.csv",
                       "05_19_2023_12_36_24_data.csv", "05_19_2023_15_45_23_data.csv",
                       "05_19_2023_10_17_29_data.csv", "05_19_2023_13_23_33_data.csv",
                       "05_19_2023_16_33_02_data.csv", "05_19_2023_11_03_17_data.csv",
                       "05_19_2023_14_10_43_data.csv", "05_19_2023_17_21_28_data.csv"
                       ]
        self.pressures = [200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 760]
        # Configure the data plots
        self.load_pressure_scan()
        
        # Default values
        sigx0, mux0, theta0, phi0 = 3.0, 6.0, 65, 0
        
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
        for ii in range(12): # FIXME (number of data files in pressure scan)
            row = int(np.floor(ii/6))
            col = ii%6

            # One set of axes per pressure
            self.dataPlotAxes.append( self.dataPlots.addPlot(row=row, col=col) )
            
            # Pressure scan data (points)
            self.dataPlotPoints.append( self.dataPlotAxes[-1].plot(self.chans[1:], self.pressure_scan_data[ii][1:],
                                                                  pen=None) )
            self.dataPlotPoints[-1].setSymbol('o')

            # "fits" to pressure scan data (lines)
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

        # total diffusion
        sigx_min, sigx_max, sigx_set = 0, 8, int(sigx0)
        self.slider_sigx = qtw.QSlider(qtc.Qt.Horizontal, self)
        self.slider_sigx.setMinimum(sigx_min)
        self.slider_sigx.setMaximum(sigx_max)
        self.slider_sigx.setValue(sigx_set)

        # Offset
        mux_min, mux_max, mux_set = 0, int(xmax), int(mux0)
        self.slider_mux = qtw.QSlider(qtc.Qt.Horizontal, self)
        self.slider_mux.setMinimum(mux_min)
        self.slider_mux.setMaximum(mux_max)
        self.slider_mux.setValue(mux_set)
        
        # Theta (angle wrt cathode normal)
        theta_min, theta_max, theta_set = 0, 85, int(theta0)
        self.slider_theta = qtw.QSlider(qtc.Qt.Horizontal, self)
        self.slider_theta.setMinimum(theta_min)
        self.slider_theta.setMaximum(theta_max)
        self.slider_theta.setValue(theta_set)

        # Phi (angle in plane of cathode)
        phi_min, phi_max, phi_set = -90, 90, int(phi0)
        self.slider_phi = qtw.QSlider(qtc.Qt.Horizontal, self)
        self.slider_phi.setMinimum(phi_min)
        self.slider_phi.setMaximum(phi_max)
        self.slider_phi.setValue(phi_set)

        # Layout
        self.layout = qtw.QVBoxLayout()

        self.layoutControls = qtw.QVBoxLayout()
        #self.controls = qtw.QWidget()
        #self.controls.setLayout(qtw.QVBoxLayout())
        self.layoutControls.addWidget(self.slider_sigx)
        self.layoutControls.addWidget(self.slider_mux)
        self.layoutControls.addWidget(self.slider_theta)
        self.layoutControls.addWidget(self.slider_phi)
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

        # Display the initial plots
        gg0 = diff.make_gaussian(self.xx, self.yy, sigx0, mux0, theta0, phi0)
        self.img.setImage(gg0.T)
        #self.img.setImage(self.masks[15])
        self.imageDisplay.setXRange(xmin, xmax)
        self.imageDisplay.setYRange(ymin, ymax)
        self.imageDisplay.setAspectLocked(True)

        # Define the transformation (from pixel number to physical units)
        tr = qtg.QTransform()  # prepare ImageItem transformation:
        tr.scale((xmax-xmin)/nx, (ymax-ymin)/ny)       # scale horizontal and vertical axes
        tr.translate( -nx/2, -ny/2 ) # locate center at axis origin
        self.img.setTransform(tr) # assign transform

        # Default profile
        self.line.setData(x=self.chans, y=self.integrateCharge(gg0))
        self.profileGraph.setTitle(f"sigx, mux, theta, phi = {sigx0}, {mux0}, {theta0}, {phi0}")

        # overlay the anode boundaries
        #p_ellipse = pg.QtGui.QGraphicsEllipseItem(0, 0, 3, 3)  # x, y, width, height
        self.draw_anode()


        self.add_profile_to_data_plots()

    def add_profile_to_data_plots(self):
        xvals, yvals = self.line.getData()
        for ii in range(len(self.fnames)):
            _, ydata = self.dataPlotPoints[ii].getData()
            # scale ydata to match actual data
            yvals *= np.max(ydata)/np.max(yvals)
            self.dataPlotLines[ii].setData(xvals[1:], yvals[1:])
       
        
    def load_pressure_scan(self):
        #dirname = "/Users/jbattat/research/qpix/saq_analysis/data/pressure_scan/20230519/"
        fnames = [os.path.join(self.dirname, fn) for fn in self.fnames]
        print(fnames)

        for ii in range(len(fnames)):
            junk1, rsts, junk2 = np.loadtxt(fnames[ii], delimiter=",", skiprows=1, unpack=True)
            self.pressure_scan_data.append(rsts[:])
        
    def draw_anode(self):
        for ii in range(16):
            w = 2*self.midpoints[ii]
            p_ellipse = qtw.QGraphicsEllipseItem(-w/2, -w/2, w, w)  # x, y, width, height
            p_ellipse.setPen(pg.mkPen((255, 0, 0, 100)))
            #p_ellipse.setBrush(pg.mkBrush((50, 50, 200)))
            self.imageDisplay.addItem(p_ellipse)
            
    def integrateCharge(self, cloud):
        ll = []
        for ii in range(16):
            ll.append( np.sum(cloud*self.masks[ii])/self.areas[ii])
        return ll
    
    def updateGraphs(self):
        sigx = self.slider_sigx.value()
        mux = self.slider_mux.value()
        theta = self.slider_theta.value()
        phi = self.slider_phi.value()
        print(f'mux, theta, phi = {mux}, {theta}, {phi}')
        gg = diff.make_gaussian(self.xx, self.yy, sigx, mux, theta, phi)
        self.img.setImage(gg.T)
        self.line.setData(x=self.chans, y=self.integrateCharge(gg))
        self.profileGraph.setTitle(f"sigx, mux, theta, phi = {sigx}, {mux}, {theta}, {phi}")

        # update lines on the pressure scan plots
        xvals, yvals = self.line.getData()
        for ii in range(len(self.fnames)):
            _, ydata = self.dataPlotPoints[ii].getData()
            # scale ydata to match actual data
            self.dataPlotLines[ii].setData(xvals[1:], yvals[1:]*np.max(ydata)/np.max(yvals))

        
if __name__ == "__main__":
    app = qtw.QApplication(sys.argv)

    win = MainWindow()
    win.show()

    sys.exit(app.exec_())
