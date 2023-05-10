# MCM2022
Code for the COMAP 2022 math contest in modeling, team control number 2226471

Team members: Trevor Loe, Madison Lytle, Callan Whitney

The competiton was completed Febuary 18-20, 2022. 

## System Specs
This code was created and run in MATLAB 2021b. It will likely work with further systems but has not been extensively tested. The optimization toolbox is required, specifically for the function ```fmincon```.

## File Description
### Simulation Scripts
- [main](main.m) This script runs the model based on the parameter values defined in [params_discrete](params_discrete.m). Saves the data to a local file in the data folder model_data/ followed by the data and time the simulation was run. This is done using the following scripts called 
  - [sqp_solve](sqp_solve.m) which makes a call to [sqp_run_new](sqp_run_new.m), the script utilizing a sequential quadratic programming method to solve the constrained optimization problem, with constraint equations given by the Skiba model. The structure ```biker``` is passed in containing all the rider parameters. The structure ```course``` is passed in, containing the necessary course parameters (the length, grade, and environmental parameters. 
