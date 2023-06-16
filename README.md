# Solar wind model

This project was given as an exercise in TFY4240 [Electromagnetic Theory](https://www.ntnu.edu/studies/courses/TFY4240#tab=omEmnet) at [NTNU](https://www.ntnu.edu/) spring 2022. 

**Project description**: Model proton trajectories in earths magnetic field

## Implementation
Lorentz force with no electric field: 

![](https://latex.codecogs.com/svg.image?{\ddot{\vec{r}}&space;=&space;\frac{q}{m}\dot{\vec{r}}\times\vec{B}})   
[Runge-Kutta-Fehlberg](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method) method was used to solve this equation.

src/main.jl : initial conditions
src/physics.jl : RKF implementation

### Example
<img src="figures/report/trajectoriesLowXZ.pdf" width="700">
