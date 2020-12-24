# BARBEQuE
**B**ring **A** **R**ational **B**ell **E**xperiment on **Qu**antum **E**xperience

A simple code that allows to verify, through the violation of **Bell's inequality** (in Wigner's version) that local realism is violated in quantum mechanics.

## Why “on Quantum Experience”?
It is possible to demonstrate the violation of Bell's inequality on a quantum computer, in this case on the [IBM Quantum Experience](https://quantum-computing.ibm.com) platform and in three different ways:
1. Local quantum simulation
2. Quantum simulation on IBM QE
3. Quantum computing on IBM QE

## The circuits
In order to write an algorithm in quantum computing, the concept of circuit is exploited. In this case, suppose we consider the entangled state of two spin 1/2 fermions: **the singlet state**

<p align="center">
  <img src="https://latex.codecogs.com/svg.latex?|00\rangle=\frac{|+-\rangle-|-+\rangle}{\sqrt2}">
</p>

In the form of a quantum circuit, exploiting the **Hadamard-gate** and **CNOT-gate** operators for the construction of this state, we have:

<p align="center">
  <img src="images/entangled.png" width="200">
</p>

In the **Wigner version** it is assumed to perform spin measurements with respect to 3 generic axes, in this case coplanars as in the following reference system where **c** bisects the angle between the **a** and **b** axes:

<p align="center">
  <img src="images/graph.png" width="400">
</p>

Carrying out the calculations it comes to the conclusion that this inequality:

<p align="center">
  <img src="https://latex.codecogs.com/svg.latex?P(+_a,+_b)\le%20P(+_a,%20+_c)+P(+_c,+_b)">
</p>

is violated for

<p align="center">
  <img src="https://latex.codecogs.com/svg.latex?\theta\in\Big(0,\frac\pi2\Big)">
</p>

The corresponding quantum circuits are as follows:

<p align="center">
  <img src="images/ab-axes.png" width="300">
</p>

Measurement of the **first qubit** along axis **a**<br>
Rotation of the **second qubit** of [![](https://latex.codecogs.com/svg.latex?2\theta)](#) and measure along axis **b**

<p align="center">
  <img src="images/ac-axes.png" width="300">
</p>

Measurement of the **first qubit** along axis **a**<br>
Rotation of the **second qubit** of [![](https://latex.codecogs.com/svg.latex?\theta)](#) and measure along axis **c**

<p align="center">
  <img src="images/cb-axes.png" width="300">
</p>

Rotation of the **first qubit** of [![](https://latex.codecogs.com/svg.latex?\theta)](#) and measure along axis **c**<br>
Rotation of the **second qubit** of [![](https://latex.codecogs.com/svg.latex?2\theta)](#) and measure along axis **b**

## Analysis
By plotting two functions as a function of the angle, it is possible to observe in which regions the inequality is satisfied (over 0) and in which it is violated (under 0).

1. Local quantum simulation

<p align="center">
  <img src="images/ls.png" width="300">
</p>

2. Quantum simulation on IBM QE

<p align="center">
  <img src="images/qs.png" width="300">
</p>

3. Quantum computing on IBM QE

<p align="center">
  <img src="images/qc.png" width="300">
</p>
