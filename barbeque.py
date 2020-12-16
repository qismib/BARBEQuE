from qiskit import *
from qiskit.visualization import *
from qiskit.tools.monitor import *
import matplotlib.pyplot as plt
import numpy as np

provider = IBMQ.load_account()

def bellExperiment(angles = [0,0], launcher = {'backend': 'qasm_simulator', 'shots': 1024}) :

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
            circ.rz(theta, index)
        circ.h(index)

    # Measuring
    circ.measure(range(2),range(2))

    # Drawing the circuit
    circ.draw(output='mpl')

    # Running the circuit on a backend

    if (launcher['backend'] == 'qasm_simulator') :
        backend = Aer.get_backend(launcher['backend'])
    else :
        backend = provider.get_backend(launcher['backend'])

    job = execute(circ, backend, shots=launcher['shots'])
    counts = job.result().get_counts(circ)

    return counts

def loadBackends() :
    # Loading the backends
    provider.backends()

    # Printing the available backends
    for b in provider.backends():
        print(b.status().to_dict())

    backend_overview()


def computeProb(launcher) :

    # Initializing
    angles = np.linspace(0, 2*np.pi, 20)
    shots = launcher['shots']
    pAB = []
    pAC = []
    pCB = []
    pAB_00 = []
    pAC_00 = []
    pCB_00 = []
    pACB_00 = []

    # Executing
    for i, angle in enumerate(angles) :
        pAB.append(bellExperiment([0, 2*angle], launcher))
        pAC.append(bellExperiment([0, angle], launcher))
        pCB.append(bellExperiment([angle, 2*angle], launcher))

        # Handling some KeyError

        try:
            pAB[i]['00']
        except KeyError:
            pAB[i]['00'] = 0

        try:
            pAB[i]['10']
        except KeyError:
            pAB[i]['10'] = 0

        try:
            pAB[i]['01']
        except KeyError:
            pAB[i]['01'] = 0

        try:
            pAB[i]['11']
        except KeyError:
            pAB[i]['11'] = 0

        try:
            pCB[i]['00']
        except KeyError:
            pCB[i]['00'] = 0

        try:
            pCB[i]['10']
        except KeyError:
            pCB[i]['10'] = 0

        try:
            pCB[i]['01']
        except KeyError:
            pCB[i]['01'] = 0

        try:
            pCB[i]['11']
        except KeyError:
            pCB[i]['11'] = 0

        try:
            pAC[i]['00']
        except KeyError:
            pAC[i]['00'] = 0

        try:
            pAC[i]['10']
        except KeyError:
            pAC[i]['10'] = 0

        try:
            pAC[i]['01']
        except KeyError:
            pAC[i]['01'] = 0

        try:
            pAC[i]['11']
        except KeyError:
            pAC[i]['11'] = 0

        # Normalizing
        pAB_00.append((pAB[i]['00'])/shots)
        pAC_00.append((pAC[i]['00'])/shots)
        pCB_00.append((pCB[i]['00'])/shots)

        # Checking Bell's inequality
        if (pAB_00[i] <= (pAC_00[i] + pCB_00[i])) :
            print("[" + str(round(angle*180/np.pi,3)) + "°]" + " La disuguaglianza è soddisfatta: " + str(round(pAB_00[i],3)) + " è <= di " + str(round(pAC_00[i] + pCB_00[i],2)))
        else :
            print("[" + str(round(angle*180/np.pi,3)) + "°]" + "La disuguaglianza NON è soddisfatta: " + str(round(pAB_00[i],3)) + " NON è <= di " + str(round(pAC_00[i] + pCB_00[i],2)))


        pACB_00.append(pAC_00[i] + pCB_00[i])

    # Plotting graphs
    plt.plot(angles*180/np.pi, pAB_00)
    plt.plot(angles*180/np.pi, pACB_00)
    plt.show()

#loadBackends()
computeProb({'backend': 'qasm_simulator', 'shots': 2048})
computeProb({'backend': 'ibmq_qasm_simulator', 'shots': 2048})
computeProb({'backend': 'ibmq_vigo', 'shots': 8192})
