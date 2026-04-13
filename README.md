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
 - run_sweep (name change to Find_optimal_constellations.mat)
 - Show beam on earth (maybe call it stationary analysis and combine it with possibilities arised with a stationary example for example UE placement?)
 - Satellite centric simulation (utilization) (somewhat example code)
    - how many UEs are connected or how large part of the area is covered by the satellite (make both plots always)



# TODO
Changes to code that i want:
- I always want to feed in a "flat array" of UE positions instead of having multiple if statements inside the coverage simulator function
- grid search needs to be more function like instead of having this weird shared control relationsship between run_sweep and grid_search. It should also be made in a way where there cannot be a non optimal solution due to randomness
- run_sweep needs to start from the heighest orbit and go down with the minimum number of satellites being exactly the previous orbit height solution. This is genius
- Equal UE distribution: I want to only have to write the region limits and number of UEs and it should just make it happen. not all this guesswork  
- optimal_constellation = load("optimal_constellations.mat"); % TODO add failsafe if not present

