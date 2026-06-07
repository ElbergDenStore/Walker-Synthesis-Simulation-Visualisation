# Constellation Synthesis Simulation and Visualisation Matlab Tool
This tool has been developed to help any constellation designer synthesise constellations, run high fidelity simulations and visualise results from day one.

## Recommended usage:
- Clone reposity
- Run numerical_walker_synthesis
  Quickly find valid constellations for your region of interest. This will make the helper function generate_default_config pull from your constellations.
- Play around with the example code
- Modify to suit your specific needs
    Modify from 100% coverage requirement to max revisit time

## Required tools
 - MATLAB tools
  - Parallel computing
  - Satellite communications
 - 

## Repository structure
- Main scripts in root
- Functions folder
  - Link budget, propagators, constellation setup, UE distribution helpers and similar
- Plotting scripts folder
  - Highly specific plotting scripts
- Validator folder
  - Precision and speed benchmarks

## Main Scripts: 
 - Numerical synthesis
 - Analytical Walker-Star (SOC)
 - run constellation simulation
    - showing how default_config, run test, plot results and show constellation wokrs

## Running on external server
For synthesis and large simulations, external servers can be used. Recommended workflow is VS code and MATLAB extension. 

matlab -nosplash -nodesktop -batch "numerical_walker_synthesis"

For some functionality a display is needed. Solution is using a virtual display via xvfb and tmux.

example: xvfb-run -a --server-args="-screen 0 1920x1080x24" matlab -nosplash -nodesktop -batch "numerical_walker_synthesis" > log.txt 

# Further work
- Remake the tx gain adjuster to use beam on earth instead of size of antenna consistent across altitudes and frequency
- Make all variables follow code style.
- Coverage simulator is ugly due to the early stopping, but it seems necessary unfortunately... I dont like it
- Add Uplink functionality
- Develop satellite focused view?
- Revive UE focus view from the dead?
