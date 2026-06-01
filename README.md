# Constellation Synthesis Simulation and Visualisation Matlab Tool
This tool has been developed to help any constellation designer to run high fidelity simulations.

Features: (I have an idea of grouping them into simulation, illustration and validation or something? hmmm i dont know yet)

It has been developed to be easily extended for future research.


## Recommended usage:
- Clone reposity
- Run Find Optimal Constellations
  Quickly find valid constellations for your region of interest. This will make the helper function generate_default_config pull from your constellations.
- Play around with the example code
- Modify to suit your specific needs
    Modify from 100% coverage requirement to max revisit time

## Required tools
 - MATLAB tools
  - Parallel computing
  - Satellite communications
 - 



## Main simulations: 
 - Numerical Walker-Delta
 - Numerical Walker-Star
 - Analytical Walker-Star (SOC)
 - run constellation simulation
    - showing how get_cfg, run test, plot results, show contellation
    - precompute
      - P618 attenuation
      - Walker-Delta constellations

Plotting and visualisation
 - Show beam on earth (maybe call it stationary analysis and combine it with possibilities arised with a stationary example for example UE placement?)
 - Interference
 - Visualise constellation

  - Trade-offs
    - antenna gain and altitude
    - 


Functions:
 - Link budget folder
    - link calc matrix
    -
 - helper functions
  - generate equal area?
 - assymmetrical walker-star construction
 - 

  beams:
  - hexagonal
  - oneweb

Validation
 - AER
 - Propagators
 - Analytical SOC and Walker-Star

## Running on external server
For synthesis and large simulations, external servers can be used. Recommended workflow is VS code and MATLAB extension. For some functionality use virtual display via xvfb and tmux.

example: xvfb-run -a --server-args="-screen 0 1920x1080x24" matlab -nosplash -nodesktop -batch "numerical_walker_synthesis" > log.txt 

# TODO
Changes to code that i want
- SEARCH FOR ALL TODOs
- Remove stupid and unnecessary comments
- Fix names of files to be more easier read
  - walker-star and walker-delta version of different files? what to do? generate_walker_delta_states.m generate_walker_star_states.m? I think it should be 1 file that has if else basically.
- Follow code style
- Big cleanup in the files
- Move all the plotting files to the plotting folder
- Coverage simulator is ugly due to the early stopping and workers might not be able to not die using it
- The gridsearch is ugly due to all the fallback mechanisms, I need to clean it up in one go.
- Create understandable read-me so others can use the code, including the dependencies.
- Develop satellite focused view?
- Add Uplink functionality
- Revive UE focus view from the dead?





option b:

main scripts in the root:
 - Numerical Walker-Delta and Walker-Star combined and combined with "Run single grid search"?
 - Analytical Walker-Star (SOC)
 - Run simulator
 - plot simulation
 - show constellation
 - get config

Plotting_scripts
Validators
Functions

Where to put generate p618 LUT? in functions?