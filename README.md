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



# FORTRAN Program to Calculate Radial Distribution Function (RDF) of Liquid Argon

## Overview

This program calculates the Radial Distribution Function for a simulated liquid Argon system. RDF is a crucial property that gives insights into the structure of a liquid by measuring how the density of particles varies as a function of distance from a reference particle. The program reads the atomic configuration from a file, computes distances between particle pairs using the minimum image convention, and calculates the RDF for the system.

## Files

1. **inputparameter.dat**: Contains parameters for the simulation, such as the number of unit cells and unit cell length.
2. **in_mcargon.dat**: Contains input data such as random seed, maximum displacement, LJ potential parameters, and temperature.
3. **finalconfi.xyz**: Holds the final configuration of atoms (positions) from a previous simulation.
4. **rdf.dat**: The output file where the calculated RDF values are stored.
5. **positionofneighbor.dat**: Logs the distances between neighboring particles and their corresponding bin values.

## Program Structure

1. **Initialization**:
   - The program reads simulation parameters from `inputparameter.dat` and `in_mcargon.dat`.
   - The size of the simulation box and the number of atoms are calculated based on the unit cell length and the number of unit cells.
   - The initial atomic configuration is read from the file `finalconfi.xyz`.

2. **Radial Distribution Function (RDF) Calculation**:
   - For each pair of atoms, the program computes the distance between them using the minimum image convention.
   - The distances are sorted into bins, each representing a specific distance interval.
   - The RDF is calculated by counting how many pairs of atoms fall within each distance bin and normalizing it by the total number of atoms and bin volume.

3. **Output**:
   - The RDF values are written to `rdf.dat` for each distance bin.
   - The distances between pairs of atoms and their corresponding bin values are written to `positionofneighbor.dat`.

## Running the Program

1. Prepare the input files (`inputparameter.dat`, `in_mcargon.dat`, `finalconfi.xyz`).
2. Compile and run the program:
   ```
   gfortran rdf_mcargon.f95 -o rdf
   ./rdf
   ```
3. The program will generate the following output files:
   - **rdf.dat**: Contains the radial distribution function for different distances.
   - **positionofneighbor.dat**: Logs the pairwise distances and the corresponding bin indices.

## Important Outputs

- **RDF (rdf.dat)**: The file contains three columns:
  - The first column is the radial distance from the reference particle.
  - The second column shows the number of particle pairs at that distance.
  - The third column gives the RDF value for that distance.
  
  The RDF describes how the probability of finding a particle at a distance `r` from a reference particle changes as `r` increases.

- **Neighbor Information (positionofneighbor.dat)**: Logs the distances between atom pairs, the corresponding bin, and normalization data for further analysis.

## Modifications

To modify the simulation:
- Adjust the bin size by modifying the divisor used in `sqrt(rsqr)/0.001d0`.
- Change the system size by modifying `nx` and `unitcelllength` in `inputparameter.dat`.
- Tweak Lennard-Jones parameters (`eps`, `sgmahex`, `sgmatwlv`) to simulate different interaction strengths between Argon atoms.

## Contact

For any issues or questions, please contact the author of the code (shubhadeepnag92@gmail.com).

