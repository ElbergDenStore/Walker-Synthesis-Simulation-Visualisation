# Satcom4Defence Constellation Simulator Matlab Tool
This tool has been developed to help any constellation designer to run high fidelity simulations.

Features: (I have an idea of grouping them into simulation, illustration and validation or something? hmmm i dont know yet)

It has been developed to be easily used for future development.


## Recommended usage:
- Clone reposity
- Run Find Optimal Constellations
  Quickly find valid constellations for your region of interest. This will make the helper function generate_default_config pull from your constellations.
- Play around with the example code
- Modify to suit your specific needs
    Modify from 100% coverage requirement to max revisit time



## Main simulations: 
 - OneWeb validation 
 - run single coverage test 
    - showing how get_cfg, run test, plot results, show contellation
 - Find_valid_constellations.mat
 - Show beam on earth (maybe call it stationary analysis and combine it with possibilities arised with a stationary example for example UE placement?)
 - Satellite centric simulation (utilization) (somewhat example code)
    - how many UEs are connected or how large part of the area is covered by the satellite (make both plots always)




# TODO
Changes to code that i want:
<!-- - OneWeb Validation should run and produce nice plots, and be correct -->
- FIX INTERFERENCE CALC
- Follow code style
- SEARCH FOR ALL TODOs
- Add Uplink functionality
- Remove stupid and unnecessary comments
- Big cleanup in the files
- Move all the plotting files to the plotting folder
- Revive UE focus view from the dead?
- Develop satellite focused view?
- I need to validate the smarter SOC with lower inclinations to also work for 90 degrees: it did so now i need to remove the old one.
- Coverage simulator is ugly due to the early stopping and workers might not be able to not die using it
- The gridsearch is ugly due to all the fallback mechanisms, I need to clean it up in one go.
- publish to github
- Create understandable read-me so others can use the code, including the dependencies.

