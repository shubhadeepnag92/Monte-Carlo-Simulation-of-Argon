# MC Program to Simulate Liquid Argon

## Overview

This program simulates the behavior of liquid Argon using the MC method. It computes the interaction energies between Argon atoms based on the Lennard-Jones (LJ) potential and uses Metropolis sampling to model the system's evolution. The simulation is carried out in a periodic simulation box with a specified number of atoms arranged in a cubic lattice. The results include the system's initial and final energies, the acceptance ratio of Monte Carlo moves, and the final configuration of atoms.

## Files

1. **inputparameter.dat**: Contains the simulation parameters such as the number of atoms and unit cell length.
2. **in_mcargon.dat**: Holds the initial seed for the random number generator and LJ potential parameters (epsilon, sigma).
3. **simu.xyz**: The input file with the initial configuration of atoms.
4. **ene.dat**: Stores the energy of the system after each cycle.
5. **output.dat**: The main output file that logs important parameters, initial and final energies, and acceptance ratios.
6. **accpprob.dat**: Stores the acceptance probability of Monte Carlo moves.
7. **finalconfi.xyz**: Contains the final configuration of atoms after the simulation.

## Simulation Parameters

- **nx**: Number of unit cells along each axis.
- **unitcelllength**: Length of each unit cell in nanometers.
- **isd**: Seed for random number generator.
- **d_m**: Maximum displacement allowed for each atom during MC moves.
- **ncyc**: Number of Monte Carlo cycles to run.
- **neq**: Number of cycles for equilibration.
- **eps**: Depth of the potential well in the Lennard-Jones potential.
- **sgmahex**: Sigma parameter for the Lennard-Jones potential (hex term).
- **sgmatwlv**: Sigma parameter for the Lennard-Jones potential (twelve term).
- **R**: Boltzmann constant.
- **T**: Temperature of the system in Kelvin.

## Program Structure

1. **Initialization**: 
   - The program reads parameters from `inputparameter.dat` and `in_mcargon.dat`.
   - It calculates the initial positions of atoms in a cubic lattice and the simulation box size.
   - The initial energy of the system is calculated based on the Lennard-Jones potential.

2. **Monte Carlo Simulation**:
   - For each Monte Carlo cycle:
     - Each atom is moved randomly within the simulation box.
     - The initial energy and new energy of the atom are computed based on its interactions with all other atoms.
     - The new configuration is accepted or rejected based on the Metropolis criterion.
   - The total energy of the system is updated after each cycle.

3. **Energy Calculations**:
   - The Lennard-Jones potential is used to compute the interaction energy between each pair of atoms.
   - Both the initial and new energy configurations are calculated, and the change in energy determines whether the move is accepted.

4. **Final Output**:
   - The final energy and atomic configuration are written to output files.
   - The acceptance probability and energy per atom are recorded.

## Running the Program

1. Prepare the necessary input files (`inputparameter.dat`, `in_mcargon.dat`, `simu.xyz`).
2. Compile and run the program:
   ```
   gfortran mcargon.f95 -o argon
   ./argon
   ```
3. The program will output the following files:
   - **output.dat**: Logs simulation details, including initial and final energies.
   - **ene.dat**: Contains the energy after each cycle.
   - **accpprob.dat**: Records acceptance probabilities.
   - **finalconfi.xyz**: Final atom positions.

## Important Outputs

- **Initial and Final Energy**: The energy of the system before and after the simulation is logged.
- **Acceptance Probability**: The ratio of accepted Monte Carlo moves to the total number of moves.
- **Final Configuration**: The atomic positions after the simulation is complete.

## Modifications

To modify the simulation:
- Change the number of cycles (`ncyc`) and maximum displacement (`d_m`) in `in_mcargon.dat`.
- Adjust the Lennard-Jones parameters (`eps`, `sgmahex`, `sgmatwlv`) to simulate different interaction strengths.
- Increase or decrease the number of atoms by modifying `nx` in `inputparameter.dat`.

## Contact

For any issues or questions, please contact the author of the code (shubhadeepnag92@gmail.com).

