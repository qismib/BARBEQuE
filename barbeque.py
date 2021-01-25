# pylint: disable=unused-wildcard-import, method-hidden
# pylint: enable=too-many-lines

import os
import json
import numpy as np
import matplotlib.pyplot as plt

from qiskit import *
from qiskit.visualization import * 
from qiskit.tools.monitor import *
from itertools import product, chain 
from scipy.optimize import least_squares
from matplotlib.ticker import FormatStrFormatter

# Fitting function with 4 params
def fitFunc(x, a, b, c):
    return a*(np.sin(b*x)**2)+a*(np.sin(b*x)**2)-c*(np.sin(x)**2)

# Residues function with 4 params
def residues(x, t, y):
    return fitFunc(t, x[0], x[1], x[2])-y

# Storing data
def storeData(name, data):

    filename = "data/%s.json" % name

    # If it exists, only update the old data
    if (os.path.exists(filename)) :
        with open(filename, "r+") as file:

            # Load the old data
            oldData = json.load(file)

            # Add the new data
            oldData.append(data)
            file.seek(0)

            # Update JSON file
            json.dump(oldData, file, indent = 4)

    # If it doesn't exist, create a JSON file
    else :
        jsonData = []
        jsonData.append(data)

        with open(filename, 'w') as file:
            json.dump(jsonData, file, indent = 4)

    return None

# Loading Data
def loadData(name):
    filename = 'data/%s.json' % name

    # If it exists, load data
    if os.path.exists(filename):
        with open(filename, 'r') as file:
            data = json.load(file)

            return data
    else:
        print('[ERROR]: Check pathname %s' % filename)
    return None

# Class to create Quantum Circuits, running them,...
class BellExperiment:

    def __init__(self, nQubit, cnot, backend, shots, mode, prob, angle):

        self.nQubit = nQubit
        self.cnot = cnot
        self.backend = backend
        self.shots = shots
        self.mode = mode
        self.prob = prob
        self.angle = angle

        self.createQC()

    def createQC(self):

        self.circ = QuantumCircuit(self.nQubit, self.nQubit)

        # Set all qubits to the state |1>
        for qubit in self.cnot:
            self.circ.x(qubit)

        # Set the entangled state
        for qubit in self.cnot:
            self.circ.h(qubit[0])
            self.circ.cx(qubit[0],qubit[1])

        # Handling the circuit
        # AB: 0 theta
        # AC: 0 2theta
        # CB: theta 2theta
        if self.mode == 0: # 2 qubits 
            if self.prob == "AB":
                self.circ.rx(2*self.angle, self.cnot[0][1])
            elif self.prob == "AC":
                self.circ.rx(self.angle, self.cnot[0][1])
            elif self.prob == "CB":
                self.circ.rx(self.angle, self.cnot[0][0])
                self.circ.rx(2*self.angle, self.cnot[0][1])
        elif self.mode == 1: # 4 qubits -> same angles
            if self.prob == "AB":
                for i in range(2):
                    self.circ.rx(2*self.angle, self.cnot[i][1])
            elif self.prob == "AC":
                for i in range(2):
                    self.circ.rx(self.angle, self.cnot[i][1])
            elif self.prob == "CB":
                for i in range(2):
                    self.circ.rx(self.angle, self.cnot[i][0])
                    self.circ.rx(2*self.angle, self.cnot[i][1])
        elif self.mode == 2: # 4 qubits -> 2nd circuit uses theta + 0.5
            if self.prob == "AB":
                self.circ.rx(2*self.angle, self.cnot[0][1])
                self.circ.rx(2*(self.angle+(0.5*np.pi/180)), self.cnot[1][1])
            elif self.prob == "AC":
                self.circ.rx(self.angle, self.cnot[0][1])
                self.circ.rx(self.angle+(0.5*np.pi/180), self.cnot[1][1])
            elif self.prob == "CB":
                self.circ.rx(self.angle, self.cnot[0][0])
                self.circ.rx(self.angle+(0.5*np.pi/180), self.cnot[1][0])
                self.circ.rx(2*self.angle, self.cnot[0][1])
                self.circ.rx(2*(self.angle+(0.5*np.pi/180)), self.cnot[1][1])

        # Measuring
        for i in range(len(self.cnot)):
            self.circ.measure(self.cnot[i], self.cnot[i]) # pylint: disable=no-member

        return None

    # Running the circuit on a backend
    def runQC(self):

        # Local simulation
        if self.backend == 'qasm_simulator':
            backend = Aer.get_backend('qasm_simulator')

        # Run on a ibmq backend
        else:
            backend = provider.get_backend(self.backend)

        # Getting data
        job = execute(self.circ, backend, shots = self.shots)
        self.counts = job.result().get_counts(self.circ)

        return None
    
    # Parsing if all index are setted and normalizing data
    def parsingData(self):

        # Get all possibile combinations of indexes
        indexNames = self.getIndexes()

        # Parsing and Normalizing
        for index in indexNames:
            if index not in self.counts:
                self.counts[index] = 0
            self.counts[index] /= self.shots

        # Create obj and name
        if len(self.cnot) == 1:

            name = 'mode_%s/%s_%sq_%s%scnot_360a_%ss_p%s' % (self.mode, self.backend, self.nQubit, self.cnot[0][0], self.cnot[0][1], self.shots, self.prob)

            if self.nQubit == 2:
                # Storing data
                storeData(name, self.counts)
            else:
                tmp = {
                    # Switch from dimension nQubit to dimension 2
                    # Remember that qubits are q4q3q2q1q0 NOT q0q1q2q3q4 (for nQubits = 5)
                    "00": sum([kel[1] for kel in self.counts.items() if self.addcts(kel[0], self.cnot[0][0], 0, self.cnot[0][1], 0)]),
                    "01": sum([kel[1] for kel in self.counts.items() if self.addcts(kel[0], self.cnot[0][0], 1, self.cnot[0][1], 0)]),
                    "10": sum([kel[1] for kel in self.counts.items() if self.addcts(kel[0], self.cnot[0][0], 0, self.cnot[0][1], 1)]),
                    "11": sum([kel[1] for kel in self.counts.items() if self.addcts(kel[0], self.cnot[0][0], 1, self.cnot[0][1], 1)])
                }
                # Storing data
                storeData(name, tmp)
        else:
            name1 = 'mode_%s/%s_%sq_%s%s%s%scnot_360a_%ss_p%s_c1' % (self.mode, self.backend, self.nQubit, self.cnot[0][0], self.cnot[0][1], self.cnot[1][0], self.cnot[1][1], self.shots, self.prob)
            name2 = 'mode_%s/%s_%sq_%s%s%s%scnot_360a_%ss_p%s_c2' % (self.mode, self.backend, self.nQubit, self.cnot[0][0], self.cnot[0][1], self.cnot[1][0], self.cnot[1][1], self.shots, self.prob)
            
            tmpc1 = {
                # Switch from dimension nQubit to dimension 2
                # Remember that qubits are q4q3q2q1q0 NOT q0q1q2q3q4 (for nQubits = 5)
                "00": sum([kel[1] for kel in self.counts.items() if self.addcts(kel[0], self.cnot[0][0], 0, self.cnot[0][1], 0)]),
                "01": sum([kel[1] for kel in self.counts.items() if self.addcts(kel[0], self.cnot[0][0], 1, self.cnot[0][1], 0)]),
                "10": sum([kel[1] for kel in self.counts.items() if self.addcts(kel[0], self.cnot[0][0], 0, self.cnot[0][1], 1)]),
                "11": sum([kel[1] for kel in self.counts.items() if self.addcts(kel[0], self.cnot[0][0], 1, self.cnot[0][1], 1)])
            }

            tmpc2 = {
                # Switch from dimension nQubit to dimension 2
                # Remember that qubits are q4q3q2q1q0 NOT q0q1q2q3q4 (for nQubits = 5)
                "00": sum([kel[1] for kel in self.counts.items() if self.addcts(kel[0], self.cnot[1][0], 0, self.cnot[1][1], 0)]),
                "01": sum([kel[1] for kel in self.counts.items() if self.addcts(kel[0], self.cnot[1][0], 1, self.cnot[1][1], 0)]),
                "10": sum([kel[1] for kel in self.counts.items() if self.addcts(kel[0], self.cnot[1][0], 0, self.cnot[1][1], 1)]),
                "11": sum([kel[1] for kel in self.counts.items() if self.addcts(kel[0], self.cnot[1][0], 1, self.cnot[1][1], 1)])
            }

            # Storing data
            storeData(name1, tmpc1)
            storeData(name2, tmpc2)

        return None

    # Getting all index combination
    def getIndexes(self):
        perm = product(['0', '1'], repeat = self.nQubit)
        possibleBin = []

        for i in list(perm):  
            myBin = ''.join(i) 
            possibleBin.append(myBin)

        return possibleBin

    # Switch from dimension nQubit to dimension 2
    def addcts(self, keystr, x, xv, y, yv):
        if (self.nQubit-1-x) > len(keystr) or (self.nQubit-1-y) > len(keystr):
            return False
        if keystr[self.nQubit-1-x] == str(xv) and keystr[self.nQubit-1-y] == str(yv):
            return True
        return False

# Class to analyzer data with plot, fit, residuals,...
class Analyzer():
    def __init__(self, nQubit, cnot, backend, shots, mode):
        self.nQubit = nQubit
        self.cnot = cnot
        self.backend = backend
        self.shots = shots
        self.mode = mode
    
    def createPlot(self):
        # Plotting figure
        self.fig = plt.figure(figsize = (9, 7.2))
        plt.frame = self.fig.add_axes((0.2, 0.4, 0.75, 0.6))

        # Labelling title, x and y
        plt.title("Bell's inequality on %s" % self.backend, fontsize = 16)
        plt.xlabel("Angles [degree]", fontsize = 12)
        plt.ylabel('Probability', fontsize = 12)

        # Labelling x-axis with step of 45 angles
        plt.xticks(np.arange(0, 361, 45))

        # Validity Limit
        plt.axhline(0, linewidth=1, linestyle="-.", color="black", label="validity-limit")
        
        # Setting the main frame
        plt.frame.get_yaxis().set_label_coords(-0.09, 0.5)
        plt.frame.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

        return None

    # Importing all data
    def importData(self):
        if self.mode == 0:

            self.pAB = loadData('mode_0/%s_%sq_%s%scnot_360a_%ss_pAB' % (self.backend, self.nQubit, self.cnot[0][0], self.cnot[0][1], self.shots))
            self.pAC = loadData('mode_0/%s_%sq_%s%scnot_360a_%ss_pAC' % (self.backend, self.nQubit, self.cnot[0][0], self.cnot[0][1], self.shots))
            self.pCB = loadData('mode_0/%s_%sq_%s%scnot_360a_%ss_pCB' % (self.backend, self.nQubit, self.cnot[0][0], self.cnot[0][1], self.shots))

        elif self.mode == 1 or self.mode == 2:

            self.pAB_c1 = loadData('mode_%s/%s_%sq_%s%s%s%scnot_360a_%ss_pAB_c1' % (self.mode, self.backend, self.nQubit, self.cnot[0][0], self.cnot[0][1], self.cnot[1][0], self.cnot[1][1], self.shots))
            self.pAC_c1 = loadData('mode_%s/%s_%sq_%s%s%s%scnot_360a_%ss_pAC_c1' % (self.mode, self.backend, self.nQubit, self.cnot[0][0], self.cnot[0][1], self.cnot[1][0], self.cnot[1][1], self.shots))
            self.pCB_c1 = loadData('mode_%s/%s_%sq_%s%s%s%scnot_360a_%ss_pCB_c1' % (self.mode, self.backend, self.nQubit, self.cnot[0][0], self.cnot[0][1], self.cnot[1][0], self.cnot[1][1], self.shots))

            self.pAB_c2 = loadData('mode_%s/%s_%sq_%s%s%s%scnot_360a_%ss_pAB_c2' % (self.mode, self.backend, self.nQubit, self.cnot[0][0], self.cnot[0][1], self.cnot[1][0], self.cnot[1][1], self.shots))
            self.pAC_c2 = loadData('mode_%s/%s_%sq_%s%s%s%scnot_360a_%ss_pAC_c2' % (self.mode, self.backend, self.nQubit, self.cnot[0][0], self.cnot[0][1], self.cnot[1][0], self.cnot[1][1], self.shots))
            self.pCB_c2 = loadData('mode_%s/%s_%sq_%s%s%s%scnot_360a_%ss_pCB_c2' % (self.mode, self.backend, self.nQubit, self.cnot[0][0], self.cnot[0][1], self.cnot[1][0], self.cnot[1][1], self.shots))

        return None

    # Evaluating a specific probability
    def evaluateProb(self, p):
        if self.mode == 0:

            self.data, pAB_00, pAC_00, pCB_00 = [[]], [], [], []

            for i in range(360):
                pAB_00.append(self.pAB[i][p])
                pAC_00.append(self.pAC[i][p])
                pCB_00.append(self.pCB[i][p])

                self.data[0].append(pAC_00[i] + pCB_00[i] - pAB_00[i])

        elif self.mode == 1 or self.mode == 2:

            self.data, pABc1, pACc1, pCBc1, pABc2, pACc2, pCBc2 = [[],[]], [], [], [], [], [], []

            for i in range(360):

                pABc1.append(self.pAB_c1[i][p])
                pACc1.append(self.pAC_c1[i][p])
                pCBc1.append(self.pCB_c1[i][p])

                self.data[0].append(pACc1[i] + pCBc1[i] - pABc1[i])
                
                pABc2.append(self.pAB_c2[i][p])
                pACc2.append(self.pAC_c2[i][p])
                pCBc2.append(self.pCB_c2[i][p])

                self.data[1].append(pACc2[i] + pCBc2[i] - pABc2[i])

        return None

    # Plotting data
    def plotData(self):
        self.angles = np.linspace(0, 2*np.pi, 360)

        if self.mode == 0:
            plt.plot(self.angles*180/np.pi, self.data[0], 'o', markersize = .7, color='#FF5600', label = "exp-data")

        elif self.mode == 1:
            plt.plot(self.angles*180/np.pi, self.data[0], 'o', markersize = .7, color='#FF5600', label = "exp-data c1")
            plt.plot(self.angles*180/np.pi, self.data[1], 'o', markersize = .7, color='#2419B2', label = "exp-data c2")

        else:
            plt.plot(self.angles*180/np.pi, self.data[0], 'o', markersize = .7, color='#FF5600', label = "exp-data c1")
            plt.plot(self.angles*180/np.pi+0.5, self.data[1], 'o', markersize = .7, color='#2419B2', label = "exp-data c2")

        return None

    # Fitting data
    def fitData(self):

        initial_guesses = np.array([0.5, 0.5, 0.5])
        
        self.fit = least_squares(residues, initial_guesses, args=(self.angles, self.data[0]), loss='cauchy', bounds=(([0,0,0], [1,1,1])))
        
        if self.mode == 0:
            plt.plot(self.angles*180/np.pi, fitFunc(self.angles, *self.fit.x), color='#FFC500', linewidth = 1, label='robust-fit')
        elif self.mode == 1:
            plt.plot(self.angles*180/np.pi, fitFunc(self.angles, *self.fit.x), color='#FFC500', linewidth = 1, label='robust-fit c1')
            
            self.fit2 = least_squares(residues, initial_guesses, args=(self.angles, self.data[1]), loss='cauchy', bounds=(([0,0,0], [1,1,1])))
            plt.plot(self.angles*180/np.pi, fitFunc(self.angles, *self.fit2.x), color='#00AA72', linewidth = 1, label='robust-fit c2')
        else:
            plt.plot(self.angles*180/np.pi, fitFunc(self.angles, *self.fit.x), color='#FFC500', linewidth = 1, label='robust-fit c1')
            
            self.fit2 = least_squares(residues, initial_guesses, args=(self.angles, self.data[1]), loss='cauchy', bounds=(([0,0,0], [1,1,1])))
            plt.plot(self.angles*180/np.pi+0.5, fitFunc(self.angles+(0.5*np.pi/180), *self.fit2.x), color='#00AA72', linewidth = 1, label='robust-fit c2')

        plt.legend(loc='upper right')

        return None

    # Plotting residuals
    def plotResiduals(self):
        plt.subframe = self.fig.add_axes((0.2,0.2,.75,.2))
        plt.xlabel("Angle [degree]", fontsize = 12)
        plt.ylabel("Residual", fontsize = 12)
        plt.xticks(np.arange(0,361,45))
        plt.subframe.get_yaxis().set_label_coords(-0.09,0.5)

        plt.axhline(0, linewidth=1, linestyle="--", color="black")

        plt.plot(self.angles*180/np.pi, self.fit.fun, "x", markersize = 1.3, color='#FF5600')

        if self.mode == 1:
            plt.plot(self.angles*180/np.pi, self.fit2.fun, "x", markersize = 1.3, color='#2419B2')

        elif self.mode == 2:
            plt.plot(self.angles*180/np.pi+0.5, self.fit2.fun, "x", markersize = 1.3, color='#2419B2')

        return None

    def showPlot(self):
        plt.show()

    # Saving plot in svg and pdf format
    def printPlot(self):
        if self.mode == 0:
            self.fig.savefig('plots/mode_%s/%s_%sq_%s%scnot_360a_%ss.svg' % (self.mode, self.backend, self.nQubit, self.cnot[0][0], self.cnot[0][1], self.shots), bbox_inches='tight')
            self.fig.savefig('plots/mode_%s/%s_%sq_%s%scnot_360a_%ss.pdf' % (self.mode, self.backend, self.nQubit, self.cnot[0][0], self.cnot[0][1], self.shots), bbox_inches='tight')
        else:
            self.fig.savefig('plots/mode_%s/%s_%sq_%s%s%s%scnot_360a_%ss.svg' % (self.mode, self.backend, self.nQubit, self.cnot[0][0], self.cnot[0][1], self.cnot[1][0], self.cnot[1][1], self.shots), bbox_inches='tight')
            self.fig.savefig('plots/mode_%s/%s_%sq_%s%s%s%scnot_360a_%ss.pdf' % (self.mode, self.backend, self.nQubit, self.cnot[0][0], self.cnot[0][1], self.cnot[1][0], self.cnot[1][1], self.shots), bbox_inches='tight')
        return None

    # Getting Residual Sum of Squares
    def getRSS(self):
        if self.mode == 0:
            return np.sum((self.fit.fun)**2)
        else:
            return {
                "c1": np.sum((self.fit.fun)**2),
                "c2": np.sum((self.fit2.fun)**2)
            }