I need anti gravity to read this file always.

I have currently made a large software tool that allows constellation engineers to simulate illustrate and optimize their constellation.
I created it to help myself to optimize a constellation, but as there are too many variables and too much complexity i cannot finish the project myself.
Therefore I would like to publish my code to github and have my master thesis product be this software tool.

This sets some requirements to the code going from my personal playground to semi professional looking code.

Software development guideline: 
Allow it to crash if the "protection code" is too complex. If I am finding all the valid indices by using >20 filter, it is possible to not have any valid indices, but only if the constellation is completely retarded or the simulation settings are fucked 

Software style guide:
I need to be consistent, I dont care about following some IEEE standard or whatever, but i want it consistent.
Cfg.Num_sats
for structs and similar, the first letter needs to be capitalized. Cfg.Orbit_height_m instead of Cfg.orbit_height_m. if not part of a struct, orbit_height_m instead of Orbit_height_m. 
PFD is an acronym for Power Flux Density. I want functions and variables to have PFD in capital letters. dBm means dB and here i want the B to be capital. Same for Gbps and MHz and similar
I should have config instead of cfg to avoid confusion. Same for other things that should not be shortene? Num Sats versus Number_Satellites? Shit i dont know what is best...
If the variable unit can be misunderstood, it needs to have a post_fix or whatever orbit_height_m, gain_dBi p_dBm and so on.
The main simulators needs to be starting with capital letter instead of the current naming
I hate the naming of AI generated code such as "%% 2. The Smart Filters". It should just be "% Filters" If I see numbering I get angry.

Plots:
The plots made by the simulator needs to be consistent as well. 
Blue means low, yellow means high
Units always presented in parantheses "Orbit Height (m)" (Capitalize all important words)
The plots should be self contained. If the simulation is using beam utilization of 20% and that changes the results by 5x, it needs to be a part of the title. If it does not change anything, the parameter should not be shown.
The plots should be saved and report ready.

Comments in code:
Docstring in all files
I dont want any stupid comments like "-FIXED fast matrix computation now", but i do want to help people going through the source code explaining what is going on.




