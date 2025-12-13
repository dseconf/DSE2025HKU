# Play leapfrogging game in class

Entry point: `leapfrog.m`, see settings in the top of the file

Main parameters: `setup.m` with many other scattered around the code

Transition probabilities for the stat-of-the-art: `trpr_matrix.m` to show on screen

Simulate paths for tech evolution: `sim_tech.m`



# Notes 

Transition probabilities are computed in `stp.m`

```
rissp=.2;

riss=1-cgrid(n)*rissp/(1+c0);  % probability that the state of the art cost does not drop 
                  % to a lower level next period

a=.2; % a parameter defining the triangular distribution for cost improvements:
      % must be in the unit interval  
```


Type of MPE the robot plays:
- `leapfrog.m`
- `game_simulation.m`