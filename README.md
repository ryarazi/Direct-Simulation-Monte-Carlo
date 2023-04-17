# Direct Simulation MonteCarlo
This is a DSMC (Direct Simulation Monte Carlo) project written in the Julia programming language.

DSMC (Direct Simulation Monte Carlo) is a computational technique used to simulate gas flows at the molecular level. It is based on the kinetic theory of gases, which describes the behavior of gases in terms of the motion of individual molecules.

In DSMC, the gas is represented as a collection of particles that move randomly and collide with each other. The simulation is divided into small time steps, during which the positions and velocities of the particles are updated and collisions are randimized by a proballistic model.

## Usage
To run the DSMC simulation, simply include the "driver.jl" file in a Julia environment.

```julia
include("driver.jl")
L = (1.e-6, 1.e-6, 1.e-6)
Ncell = (20, 1, 1)
ρ = 1.78
Nsim = 2000
v_init = 400.
sys = System(L, Ncell, gas, ρ, Nsim, v_init)
samp = simulate_dsmc(sys=sys, Nsteps=2000, Nsample=100)
```

This will run the simulation and return a Samle object which is a combination of samples of the system at different times.

## Configuration
The simulation can be configured the input to simulate_dsmc. This function contains the following parameters:

$L$: A $D$-dimensions tuple with the length of the system at each dimension.

$Ncell$: A $D$-dimensions tuple with the number of cells in each dimension.

gas: An object which contains the physical data of the simulated gas (see "gas.jl).

$\rho$: The Density of the System.

$Nsim$: Number of simulated particles.

v_init: The initial velocity of each particle in every axis.

$Nsteps$: Number of steps in the simulation.

$Nsample$: The sample rate for the system.

## License
This project is licensed under the MIT License. See the "LICENSE" file for details.