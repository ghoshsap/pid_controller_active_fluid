# Propotional Integral controller to switch average speed of active fluid sample on demand

PDE solver and dependencies installation: [cuPSS](https://github.com/fcaballerop/cuPSS). <br>

## PID controller setup
<img src="https://github.com/ghoshsap/pid_controller_active_fluid/blob/main/init/method.png" alt="Diagram" width="400" />

## Description 

- `pid_solver.cpp`: This file initiates integration and control process starting from an aligned initial condition.

- `pid_solver_1.cpp`: This file initiates integration and control process from arbitrary initial conditions, allowing for more varied starting points in simulations. Two example initial conditions are provided in the `init` folder.
  
- `wrapper.sh`: make it into a executable to submit a parameter sweep of jobs using `run_sim.sh`.  
