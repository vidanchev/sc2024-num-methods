# SC2024 Bootcamp Numerical Methods Workshop

## Introduction

This is some reference code for numerical methods to be used for Space Challenges 2024 Bootcamp workshops in Aug 7-8 2024. The purpose of the workshop is to familiarize cadets with numerical methods for solving differential equations and to write a basic orbital propagator which can be compared to tools such as General Mission Analysis Tool (GMAT) and Systems Toolkit (STK). The workshop is meant to be conducted as 2 sessions per ~2-3 hours with the following subdivision of tasks during each hour:

- Session 1: Intro to differential equations:
  - what do DEs describe and how do they look
  - what is the difference between continuus and discrete systems and how do we solve DEs with computers the simplest concept of solving ODE: Euler method
- Session 2: Hands-on solving of an oscillator:
  - writing an Euler integrator for a simple harmonic oscillator problem
  - plotting the Euler integrator solution vs. the real solution
  - extending to an intermediate term - the Verlet Method, comparing results with Euler (HW)
- Session 3: Extending 1D oscillator to 2D projectile motion
  - going full throttle - RK(4) method, comparing with Verlet and Euler
  - introducing the more complex right-hand side for orbital propagator or ballistic trajectory - dimensionless quantities
- Session 4: Comparing the results of our propagator with the GMAT/STK one
  - handling multiple dimensions in the RHS function
  - writing a Verlet and RK(4) method integration for 2D ballistic trajectory (with drag)
  - making your own 3D orbital propagator - HW to compare with STK data

## Files Included

- README.md the readme file
- main.py includes the main driver for the code where we set up the initial conditions and actually run the integrator which will be outputting in .csv files.
- Integrators.py includes the numerical integrators and differential equations right hand sides as functions which will be called from the main code.
- Plotter.py has a parser for reading the results and different plotting functions for the various files we're going to be exporting from the solvers.
- Satellite_PVT_GMAT.csv is the .csv file with solution after propagating the trajectory for 1 day with GMAT/STK with higher order integrator and higher J terms in the potential - we wll be comparing to this file
- Satellite_PVT_WS.csv is the .csv file with the solution after propagating the trajectory for 1 day with the RK4 propagator and only Newtonian term