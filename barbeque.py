###########################################################
# BARBEQuE                                                #
# Bring A Rational Bell Experiment on Quantum Experience  #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Repository: https://github.com/qismib/BARBEQuE          #
# Marco Gobbo: https://github.com/marcogobbo              #
# Organization: https://github.com/qismib                 #
###########################################################

import os
import json
import numpy as np

from qiskit import *
from qiskit.visualization import *
from qiskit.tools.monitor import *

# Debug Mode
DEBUG_MODE = False

# Set IBMQ API Token
IBMQ.save_account('YOUR_API_KEY', overwrite=True)

provider = IBMQ.load_account()

# Fitting function with 4 params
def fitFunc(x, a, b, c):
    return a*(np.sin(b*x)**2)+a*(np.sin(b*x)**2)-c*(np.sin(x)**2)

# Load all backends
def loadBackends():

    provider.backends()

    for b in provider.backends():
        print(b.status().to_dict())

    backend_overview()

    if DEBUG_MODE:
        print("[INFO]: Backends have been loaded.")

    return None

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

    return None

# Read data saved
def loadData(name):
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

    return None

# Bell's Experiment that create quantum circuit, run it on a backends and store the data
def bellExperiment(name, angles, launcher, mode):

    # Running the circuit on a backend
    def runCircuit():

        if launcher['backend'] == 'qasm_simulator':
            backend = Aer.get_backend('qasm_simulator')
        else:
            backend = provider.get_backend(launcher['backend'])

        job = execute(circ, backend, shots = launcher['shots'])
        counts = job.result().get_counts(circ)

        return counts

    # Mode 0: 2 qubits
    if mode == 0:

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

        # Run it
        counts = runCircuit()

        if DEBUG_MODE:
            print("[INFO]: Quantum Circuit has been executed.")

        # Parsing if all index are setted and normalizing data
        indexNames = ["00", "01", "10", "11"]

        for index in indexNames:
            if index not in counts:
                counts[index] = 0
                if DEBUG_MODE:
                    print("[INFO]: %s has been created." % index)

            # Normalizing
            counts[index] = counts[index]/launcher['shots']

        # Storing data
        storeData("mode_0/%s" % name, counts)

    # Mode 1-2: 4 qubits in parallel (same angle for each entangled state or half angles for each entangled state)
    elif mode == 1 or mode == 2:

        circ = QuantumCircuit(4,4)

        # Starting from the state |1> for all qubits
        circ.x(0)
        circ.x(1)
        circ.x(2)
        circ.x(3)

        # Creating entangled states
        circ.h(0)
        circ.cx(0,1)
        circ.h(2)
        circ.cx(2,3)

        # Doing a few rotations if required
        # Mode 1 => Same angles
        if mode == 1:
            for index, theta in enumerate(angles):
                if (theta != 0) :
                    circ.rx(theta, index)
                    circ.rx(theta, index+2)
        # Mode 2 => Half angles for each entangled state
        else:
            print("[WORK IN PROGRESS]: Need to find a solution to split angles")

        # Measuring
        circ.measure(range(4),range(4))

        # Run it
        counts = runCircuit()

        if DEBUG_MODE:
            print("[INFO]: Quantum Circuit has been executed.")

        # Parsing if all index are setted and normalizing data
        indexNames = ["0000", "0001", "0010", "0011", "0100", "0101", "0110", "0111", "1000", "1001", "1010", "1011", "1100", "1101", "1110", "1111"]

        for index in indexNames:
            if index not in counts:
                counts[index] = 0
                if DEBUG_MODE:
                    print("[INFO]: %s has been created." % index)

            # Normalizing
            counts[index] = counts[index]/launcher['shots']

        countsC1 = {
            "00": counts["0000"] + counts["0100"] + counts["1000"] + counts["1100"],
            "01": counts["0001"] + counts["0101"] + counts["1001"] + counts["1101"],
            "10": counts["0010"] + counts["0110"] + counts["1010"] + counts["1110"],
            "11": counts["0011"] + counts["0111"] + counts["1011"] + counts["1111"]
        }

        countsC2 = {
            "00": counts["0000"] + counts["0001"] + counts["0010"] + counts["0011"],
            "01": counts["0100"] + counts["0101"] + counts["0110"] + counts["0111"],
            "10": counts["1000"] + counts["1001"] + counts["1010"] + counts["1011"],
            "11": counts["1100"] + counts["1101"] + counts["1110"] + counts["1111"]
        }

        # Storing data
        if mode == 1:
            storeData("mode_1/%s_c1" % name, countsC1)
            storeData("mode_1/%s_c2" % name, countsC2)
        else:
            storeData("mode_2/%s_c1" % name, countsC1)
            storeData("mode_2/%s_c2" % name, countsC2)

    else:
        print("[ERROR]: 'mode' can assume only the values 0, 1 and 2")

    return None
