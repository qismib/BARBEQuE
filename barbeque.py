# pyright: reportUndefinedVariable=false

import os
import json
import numpy as np
import matplotlib.pyplot as plt

from qiskit import *
from qiskit.visualization import * 
from qiskit.tools.monitor import *
from itertools import product
from scipy.optimize import least_squares
from matplotlib.ticker import FormatStrFormatter

# Fitting function with 4 params
def fit_func(x, a, b, c):
    return a*(np.sin(b*x)**2)+a*(np.sin(b*x)**2)-c*(np.sin(x)**2)

# Residues function with 4 params
def residues(x, t, y):
    return fit_func(t, x[0], x[1], x[2])-y

# Storing data
def store_data(name, data):

    filename = "data/%s.json" % name

    # If it exists, only update the old data
    if (os.path.exists(filename)) :
        with open(filename, "r+") as file:

            # Load the old data
            old_data = json.load(file)

            # Add the new data
            old_data.append(data)
            file.seek(0)

            # Update JSON file
            json.dump(old_data, file, indent = 4)

    # If it doesn't exist, create a JSON file
    else :
        json_data = []
        json_data.append(data)

        with open(filename, 'w') as file:
            json.dump(json_data, file, indent = 4)

    return None

# Loading Data
def load_data(name):
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

    def __init__(self, n_qubit, cnot, backend, shots, mode, prob, angle):

        self.n_qubit = n_qubit
        self.cnot = cnot
        self.backend = backend
        self.shots = shots
        self.mode = mode
        self.prob = prob
        self.angle = angle

        self.create_qc()

    def create_qc(self):

        self.circ = QuantumCircuit(self.n_qubit, self.n_qubit)

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
            self.circ.measure(self.cnot[i], self.cnot[i])

        return None

    # Running the circuit on a backend
    def run_qc(self):

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
    def parsing_data(self):

        # Get all possibile combinations of indexes
        index_names = self.get_indexes()

        # Parsing and Normalizing
        for index in index_names:
            if index not in self.counts:
                self.counts[index] = 0
            self.counts[index] /= self.shots

        # Create obj and name
        if len(self.cnot) == 1:

            name = 'mode_%s/%s_%sq_%s%scnot_360a_%ss_p%s' % (self.mode, self.backend, self.n_qubit, self.cnot[0][0], self.cnot[0][1], self.shots, self.prob)

            if self.n_qubit == 2:
                # Storing data
                store_data(name, self.counts)
            else:
                tmp = {
                    # Switch from dimension n_qubit to dimension 2
                    # Remember that qubits are q4q3q2q1q0 NOT q0q1q2q3q4 (for n_qubits = 5)
                    "00": sum([kel[1] for kel in self.counts.items() if self.addcts(kel[0], self.cnot[0][0], 0, self.cnot[0][1], 0)]),
                    "01": sum([kel[1] for kel in self.counts.items() if self.addcts(kel[0], self.cnot[0][0], 1, self.cnot[0][1], 0)]),
                    "10": sum([kel[1] for kel in self.counts.items() if self.addcts(kel[0], self.cnot[0][0], 0, self.cnot[0][1], 1)]),
                    "11": sum([kel[1] for kel in self.counts.items() if self.addcts(kel[0], self.cnot[0][0], 1, self.cnot[0][1], 1)])
                }
                # Storing data
                store_data(name, tmp)
        else:
            name1 = 'mode_%s/%s_%sq_%s%s%s%scnot_360a_%ss_p%s_c1' % (self.mode, self.backend, self.n_qubit, self.cnot[0][0], self.cnot[0][1], self.cnot[1][0], self.cnot[1][1], self.shots, self.prob)
            name2 = 'mode_%s/%s_%sq_%s%s%s%scnot_360a_%ss_p%s_c2' % (self.mode, self.backend, self.n_qubit, self.cnot[0][0], self.cnot[0][1], self.cnot[1][0], self.cnot[1][1], self.shots, self.prob)
            
            tmpc1 = {
                # Switch from dimension n_qubit to dimension 2
                # Remember that qubits are q4q3q2q1q0 NOT q0q1q2q3q4 (for n_qubits = 5)
                "00": sum([kel[1] for kel in self.counts.items() if self.addcts(kel[0], self.cnot[0][0], 0, self.cnot[0][1], 0)]),
                "01": sum([kel[1] for kel in self.counts.items() if self.addcts(kel[0], self.cnot[0][0], 1, self.cnot[0][1], 0)]),
                "10": sum([kel[1] for kel in self.counts.items() if self.addcts(kel[0], self.cnot[0][0], 0, self.cnot[0][1], 1)]),
                "11": sum([kel[1] for kel in self.counts.items() if self.addcts(kel[0], self.cnot[0][0], 1, self.cnot[0][1], 1)])
            }

            tmpc2 = {
                # Switch from dimension n_qubit to dimension 2
                # Remember that qubits are q4q3q2q1q0 NOT q0q1q2q3q4 (for n_qubits = 5)
                "00": sum([kel[1] for kel in self.counts.items() if self.addcts(kel[0], self.cnot[1][0], 0, self.cnot[1][1], 0)]),
                "01": sum([kel[1] for kel in self.counts.items() if self.addcts(kel[0], self.cnot[1][0], 1, self.cnot[1][1], 0)]),
                "10": sum([kel[1] for kel in self.counts.items() if self.addcts(kel[0], self.cnot[1][0], 0, self.cnot[1][1], 1)]),
                "11": sum([kel[1] for kel in self.counts.items() if self.addcts(kel[0], self.cnot[1][0], 1, self.cnot[1][1], 1)])
            }

            # Storing data
            store_data(name1, tmpc1)
            store_data(name2, tmpc2)

        return None

    # Getting all index combination
    def get_indexes(self):
        perm = product(['0', '1'], repeat = self.n_qubit)
        possible_bin = []

        for i in list(perm):  
            my_bin = ''.join(i) 
            possible_bin.append(my_bin)

        return possible_bin

    # Switch from dimension n_qubit to dimension 2
    def addcts(self, keystr, x, xv, y, yv):
        if (self.n_qubit-1-x) > len(keystr) or (self.n_qubit-1-y) > len(keystr):
            return False
        if keystr[self.n_qubit-1-x] == str(xv) and keystr[self.n_qubit-1-y] == str(yv):
            return True
        return False

# Class to analyzer data with plot, fit, residuals,...
class Analyzer():
    def __init__(self, n_qubit, cnot, backend, shots, mode):
        self.n_qubit = n_qubit
        self.cnot = cnot
        self.backend = backend
        self.shots = shots
        self.mode = mode
    
    def create_plot(self):
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
        plt.axvspan(0, 90, color='lightgrey', alpha=0.5, lw=0)
        plt.axvspan(270, 360, color='lightgrey', alpha=0.5, lw=0)

        if self.mode == 0:
            angles = np.linspace(0, 360, 360)

            data_array = np.array(self.data[0])
            pos1 = np.where(data_array[:180]<0)[0]
            pos2 = np.where(data_array[180:]<0)[0] + 180

            plt.axvline(angles[pos1][0], color='#FF5600', ls='--')
            plt.axvline(angles[pos1][-1], color='#FF5600', ls='--')

            plt.axvline(angles[pos2][0], color='#FF5600', ls='--')
            plt.axvline(angles[pos2][-1], color='#FF5600', ls='--')
        
        else:
            angles = np.linspace(0, 360, 360)

            data_array1 = np.array(self.data[0])
            pos1 = np.where(data_array1[:180]<0)[0]
            pos2 = np.where(data_array1[180:]<0)[0] + 180

            try:
                plt.axvline(angles[pos1][0], color='#FF5600', ls='--')
                plt.axvline(angles[pos1][-1], color='#FF5600', ls='--')
            except:
                pass
            
            try:
                plt.axvline(angles[pos2][0], color='#FF5600', ls='--')
                plt.axvline(angles[pos2][-1], color='#FF5600', ls='--')
            except:
                pass

            data_array2 = np.array(self.data[1])
            pos3 = np.where(data_array2[:180]<0)[0]
            pos4 = np.where(data_array2[180:]<0)[0] + 180

            try:
                plt.axvline(angles[pos3][0], color='#2419B2', ls='--')
                plt.axvline(angles[pos3][-1], color='#2419B2', ls='--')
            except:
                pass

            try:
                plt.axvline(angles[pos4][0], color='#2419B2', ls='--')
                plt.axvline(angles[pos4][-1], color='#2419B2', ls='--')
            except:
                pass

        # Setting the main frame
        plt.frame.get_yaxis().set_label_coords(-0.09, 0.5)
        plt.frame.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

        return None

    # Importing all data
    def import_data(self):
        if self.mode == 0:

            self.pAB = load_data('mode_0/%s_%sq_%s%scnot_360a_%ss_pAB' % (self.backend, self.n_qubit, self.cnot[0][0], self.cnot[0][1], self.shots))
            self.pAC = load_data('mode_0/%s_%sq_%s%scnot_360a_%ss_pAC' % (self.backend, self.n_qubit, self.cnot[0][0], self.cnot[0][1], self.shots))
            self.pCB = load_data('mode_0/%s_%sq_%s%scnot_360a_%ss_pCB' % (self.backend, self.n_qubit, self.cnot[0][0], self.cnot[0][1], self.shots))

        elif self.mode == 1 or self.mode == 2:

            self.pAB_c1 = load_data('mode_%s/%s_%sq_%s%s%s%scnot_360a_%ss_pAB_c1' % (self.mode, self.backend, self.n_qubit, self.cnot[0][0], self.cnot[0][1], self.cnot[1][0], self.cnot[1][1], self.shots))
            self.pAC_c1 = load_data('mode_%s/%s_%sq_%s%s%s%scnot_360a_%ss_pAC_c1' % (self.mode, self.backend, self.n_qubit, self.cnot[0][0], self.cnot[0][1], self.cnot[1][0], self.cnot[1][1], self.shots))
            self.pCB_c1 = load_data('mode_%s/%s_%sq_%s%s%s%scnot_360a_%ss_pCB_c1' % (self.mode, self.backend, self.n_qubit, self.cnot[0][0], self.cnot[0][1], self.cnot[1][0], self.cnot[1][1], self.shots))

            self.pAB_c2 = load_data('mode_%s/%s_%sq_%s%s%s%scnot_360a_%ss_pAB_c2' % (self.mode, self.backend, self.n_qubit, self.cnot[0][0], self.cnot[0][1], self.cnot[1][0], self.cnot[1][1], self.shots))
            self.pAC_c2 = load_data('mode_%s/%s_%sq_%s%s%s%scnot_360a_%ss_pAC_c2' % (self.mode, self.backend, self.n_qubit, self.cnot[0][0], self.cnot[0][1], self.cnot[1][0], self.cnot[1][1], self.shots))
            self.pCB_c2 = load_data('mode_%s/%s_%sq_%s%s%s%scnot_360a_%ss_pCB_c2' % (self.mode, self.backend, self.n_qubit, self.cnot[0][0], self.cnot[0][1], self.cnot[1][0], self.cnot[1][1], self.shots))

        return None

    # Evaluating a specific probability
    def evaluate_prob(self, p):
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
    def plot_data(self):
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
    def fit_data(self):

        initial_guesses = np.array([0.5, 0.5, 0.5])
        
        self.fit = least_squares(residues, initial_guesses, args=(self.angles, self.data[0]), loss='cauchy', bounds=(([0,0,0], [1,1,1])))
        
        if self.mode == 0:
            plt.plot(self.angles*180/np.pi, fit_func(self.angles, *self.fit.x), color='#FFC500', linewidth = 1, label='robust-fit')
        elif self.mode == 1:
            plt.plot(self.angles*180/np.pi, fit_func(self.angles, *self.fit.x), color='#FFC500', linewidth = 1, label='robust-fit c1')
            
            self.fit2 = least_squares(residues, initial_guesses, args=(self.angles, self.data[1]), loss='cauchy', bounds=(([0,0,0], [1,1,1])))
            plt.plot(self.angles*180/np.pi, fit_func(self.angles, *self.fit2.x), color='#00AA72', linewidth = 1, label='robust-fit c2')
        else:
            plt.plot(self.angles*180/np.pi, fit_func(self.angles, *self.fit.x), color='#FFC500', linewidth = 1, label='robust-fit c1')
            
            self.fit2 = least_squares(residues, initial_guesses, args=(self.angles, self.data[1]), loss='cauchy', bounds=(([0,0,0], [1,1,1])))
            plt.plot(self.angles*180/np.pi+0.5, fit_func(self.angles+(0.5*np.pi/180), *self.fit2.x), color='#00AA72', linewidth = 1, label='robust-fit c2')

        plt.legend(loc='upper right')

        return None

    # Plotting residuals
    def plot_residuals(self):
        plt.subframe = self.fig.add_axes((0.2,0.2,.75,.2))
        plt.xlabel("Angle [degree]", fontsize = 12)
        plt.ylabel("Residual", fontsize = 12)
        plt.xticks(np.arange(0,361,45))
        plt.subframe.get_yaxis().set_label_coords(-0.09,0.5)

        plt.plot(self.angles*180/np.pi, self.fit.fun, "x", markersize = 1.3, color='#FF5600')

        if self.mode == 1:
            plt.plot(self.angles*180/np.pi, self.fit2.fun, "x", markersize = 1.3, color='#2419B2')

        elif self.mode == 2:
            plt.plot(self.angles*180/np.pi+0.5, self.fit2.fun, "x", markersize = 1.3, color='#2419B2')

        return None

    def show_plot(self):
        plt.show()

    # Saving plot in svg and pdf format
    def print_plot(self):
        if self.mode == 0:
            self.fig.savefig('plots/mode_%s/%s_%sq_%s%scnot_360a_%ss.svg' % (self.mode, self.backend, self.n_qubit, self.cnot[0][0], self.cnot[0][1], self.shots), bbox_inches='tight')
            self.fig.savefig('plots/mode_%s/%s_%sq_%s%scnot_360a_%ss.pdf' % (self.mode, self.backend, self.n_qubit, self.cnot[0][0], self.cnot[0][1], self.shots), bbox_inches='tight')
        else:
            self.fig.savefig('plots/mode_%s/%s_%sq_%s%s%s%scnot_360a_%ss.svg' % (self.mode, self.backend, self.n_qubit, self.cnot[0][0], self.cnot[0][1], self.cnot[1][0], self.cnot[1][1], self.shots), bbox_inches='tight')
            self.fig.savefig('plots/mode_%s/%s_%sq_%s%s%s%scnot_360a_%ss.pdf' % (self.mode, self.backend, self.n_qubit, self.cnot[0][0], self.cnot[0][1], self.cnot[1][0], self.cnot[1][1], self.shots), bbox_inches='tight')
        return None

    # Getting Residual Sum of Squares
    def get_rss(self):
        if self.mode == 0:
            return np.sum((self.fit.fun)**2)
        else:
            return {
                "c1": np.sum((self.fit.fun)**2),
                "c2": np.sum((self.fit2.fun)**2)
            }