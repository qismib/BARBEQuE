###########################################################
# BARBEQuE                                                #
# Bring A Rational Bell Experiment on Quantum Experience  #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Repository: https://github.com/qismib/BARBEQuE          #
# Marco Gobbo: https://github.com/marcogobbo              #
# Organization: https://github.com/qismib                 #
###########################################################

%matplotlib inline

import os
import json
import numpy as np
import matplotlib.pyplot as plt

from qiskit import *
from qiskit.visualization import *
from qiskit.tools.monitor import *
from scipy.optimize import curve_fit

# Debug Mode
DEBUG_MODE = False

# Set IBMQ API Token
provider = IBMQ.enable_account('MY_API_TOKEN')

# Load all backends
def loadBackends():

    provider.backends()

    for b in provider.backends():
        print(b.status().to_dict())

    backend_overview()

    if DEBUG_MODE:
        print("[INFO]: Backends have been loaded.")

# Store all data in JSON format
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

        if DEBUG_MODE:
            print("[INFO]: %s.json has been updated." % name)

    # If it doesn't exist, create a JSON file
    else :
        jsonData = []
        jsonData.append(data)

        with open(filename, 'w') as file:
            json.dump(jsonData, file, indent = 4)

        if DEBUG_MODE:
            print("[INFO]: %s.json has been created." % name)

# Read data saved
def readingData(name):

    filename = "data/%s.json" % name

    if os.path.exists(filename):
        with open(filename, "r") as file:
            data = json.load(file)

            if DEBUG_MODE:
                print("[INFO]: Data from %s.json has been loaded." %name)

            return data
    else :
        if DEBUG_MODE:
            print("[ERROR]: Filename %s.json doesn't exist." % name)


# Bell's Experiment that create a quantum circuit, run it on a backends and store the data
def bellExperiment(name, angles, launcher):

    # Define 2 quantum and 2 classical bit
    circ = QuantumCircuit(2,2)

    # Starting from the state |1> for both qubits
    circ.x(0)
    circ.x(1)

    # Creating an entangled state
    circ.h(0)
    circ.cx(0,1)

    # Doing a few rotations if required
    for index, theta in enumerate(angles):
        if (theta != 0) :
            circ.rx(theta, index)

    # Measuring
    circ.measure(range(2),range(2))

    # Running the circuit on a backend
    if (launcher['backend'] == 'qasm_simulator') :
        backend = Aer.get_backend('qasm_simulator')
    else :
        backend = provider.get_backend(launcher['backend'])

    job = execute(circ, backend, shots=launcher['shots'])
    counts = job.result().get_counts(circ)

    if DEBUG_MODE:
        print("[INFO]: Quantum Circuit has been executed.")

    # Parsing if all index are setted and normalizing data
    indexNames = ["00", "01", "10", "11"]

    for index in indexNames:
        if index not in counts :
            counts[index] = 0
            if DEBUG_MODE:
                print("[INFO]: %s has been created." % index)

        # Normalizing
        counts[index] = counts[index]/launcher['shots']

    # Storing data
    storeData(name, counts)

# Compute probabilities pAB, pAC, pCB
def computeProb(settings) :

    angles = settings['angles']
    backend = settings['launcher']['backend']
    shots = settings['launcher']['shots']

    # Executing
    for i, angle in enumerate(angles) :
        bellExperiment("%s_%sd%ss_pAB" % (backend, len(angles), shots), [0, 2*angle], settings['launcher'])
        bellExperiment("%s_%sd%ss_pAC" % (backend, len(angles), shots), [0, angle], settings['launcher'])
        bellExperiment("%s_%sd%ss_pCB" % (backend, len(angles), shots), [angle, 2*angle], settings['launcher'])

    if DEBUG_MODE:
        print("[INFO]: The compute has been completed successfully.")

# Fitting function with 4 params
def fitFunc(x, a, b, c) :
    return a*(np.sin(b*x)**2)+a*(np.sin(b*x)**2)-c*(np.sin(x)**2)

# Analyze the data
def analyzer(settings) :

    data, pAB_00, pAC_00, pCB_00 = [], [], [], []

    angles = settings['angles']
    backend = settings['launcher']['backend']
    shots = settings['launcher']['shots']

    # Reading data
    pAB = readingData("%s_%sd%ss_pAB" % (backend, len(angles), shots))
    pAC = readingData("%s_%sd%ss_pAC" % (backend, len(angles), shots))
    pCB = readingData("%s_%sd%ss_pCB" % (backend, len(angles), shots))

    if DEBUG_MODE:
        print("[INFO]: All data has been loaded.")

    # Evaluating only 00
    for i, angle in enumerate(angles) :
        pAB_00.append(pAB[i]['00'])
        pAC_00.append(pAC[i]['00'])
        pCB_00.append(pCB[i]['00'])

        data.append(pAC_00[i] + pCB_00[i] - pAB_00[i])

    if DEBUG_MODE:
        print("[INFO]: P(00) has been evaluated.")

    # Plotting figure
    fig = plt.figure()

    # Const function (Classical limit)
    x = np.arange(angles[-1]*180/np.pi)
    plt.plot(x, x*0)

    # Data distribution (Data processed)
    plt.plot(angles*180/np.pi, data, 'ro', markersize = 1, label = "experimental-data")

    if DEBUG_MODE:
        print("[INFO]: All data has been plotted.")

    # Initial guess (Theoretical value)
    initial_guess = [0.5, 0.5, 0.5]

    # Fitting data
    popt, pcov = curve_fit(fitFunc, angles, data, initial_guess)
    plt.plot(angles*180/np.pi, fitFunc(angles, *popt), 'y', label='fit-params: a=%5.2f, b=%5.2f, c=%5.2f' % tuple(popt))

    if DEBUG_MODE:
        print("[INFO]: Fit has been done.")

    # Labelling axes
    plt.xlabel('Angle [degree]')
    plt.ylabel('Probability')

    # Plotting graph
    plt.title("Bell's inequality on %s" % backend)
    plt.legend(loc='best')
    plt.show()

    # Print the plot
    fig.savefig("plots/%s_%sd%ss.pdf" % (backend, len(angles), shots) , bbox_inches='tight')

    if DEBUG_MODE:
        print("[INFO]: Plot has been printed.")
        print("[INFO]: Everything works. ")

# loadBackends()
# computeProb({'angles': np.linspace(0, 2*np.pi, 100), 'launcher': {'backend': 'qasm_simulator', 'shots': 10000}})
# analyzer({'angles': np.linspace(0, 2*np.pi, 100), 'launcher': {'backend': 'qasm_simulator', 'shots': 10000}})
