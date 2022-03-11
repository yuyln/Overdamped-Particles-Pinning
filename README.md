# Overdamped Particles Simulation With Pinnings Arrays

This code simulates overdamped particles that follows the Langevin equation:

$$
\alpha_d \mathbf{v}_i+\alpha_m \hat{\mathbf{z}}\times\mathbf{v}_i=\sum_{i\neq j}\mathbf{F}^{PP}_{ij}+\sum_{j}\mathbf{F}^{Pp}_{ij}+\mathbf{F}^{C}
$$

Where $\mathbf{v}_i$ is the particle velocity, $\alpha_d$ is the damping constant, that arises from
dissipation processes. $\alpha_m$ is the Magnus term, very important for some overdamped particles, such as magnetic skyrmions.

The right side of the equation are forces: $\mathbf{F}^{PP}_{ij}$ is the particle-particle interaction,
$\mathbf{F}^{Pp}_{ij}$ is the particle-pinning interaction and $\mathbf{F}^{C}$ is the driven currents force.

For the particle-particle, I modeled the potential as $U^{PP}_{ij}=U_jK_0(r_{ij})$, where $K_0$ is the modified bessel of first kind, $U_j$ is the potential strength of particle $j$, and $r_{ij}$ is the distance between particles $i$ and $j$. From this I get the force $\mathbf{F}^{PP}_{ij}=-\boldsymbol{\nabla}U^{PP}_{ij}=U_jK_1(r_{ij})\hat{\mathbf{r}}_{ij}$ where $K_1$ is the modified bessel of second kind.

For the particle-pinning I choose a gaussian potential, with the form $U^{Pp}_{ij}=U_0\exp\left(\dfrac{r_{ij}^2}{a_0^2}\right)$ where $r_{ij}$ is the distance between particle $i$ and pinning $j$, $U_0$ is the pinning strength and $a_0$ is the pinning radius. From the potential I get the force $\mathbf{F}^{Pp}_{ij}=-\boldsymbol{\nabla}U^{Pp}_{ij}=F_0r_{ij}\exp\left(\dfrac{r_{ij}^2}{a_0^2}\right)\hat{\mathbf{r}}_{ij}$ where $F_0=\dfrac{2U_0}{a_0^2}$.

The last term, is the current force, which can be of any form, in special here I consider the most general form of $\mathbf{F}^{C}=\mathbf{F}^{DC}+\mathbf{F}^{AC}$. For more detail see the input file `input/input.in`, where all the options are avaiable.

I used many molecular dynamics technics for optimization, such as subboxes, mirror boxes, lookup tables. However, if you have any further optimization, please send my a [email](mailto:jc.souza@unesp.br)

Any overdamped particle can be simulated using this code, with minor or no, changes.

This code uses Simplified Generalized Simulated Annealing for finding the ground configuration.

You can choose any number of particles and pinnings, just putting the files on `input/particles`
or `input/pinnings`.

The particles file has the following template:

```
Beta over Damping BETADAMP: VALUE
Interaction Force U0: VALUE
Positions:
10.233721388678264	18.495167551880364
.                   .
.                   .
.                   .
```
Here, `U0` is the potential strength, higher values means the particles will be more repulsive, and
lower values less repulsive (if `U0<0` then the particles will be attractive). `BETADAMP` is the parameter from the Magnus force. This is very useful for Skyrmions, for example, if your overdamped particle doesn't have Magnus force, just put `0.0`. The `positions` is the 2D positions for your particles. If you want particles of different species, just put another file with the new species data, and the program will simulate both interacting. One more thing about `BETADAMP`, the normalization $\color{lightgray}\alpha_m^2+\alpha_d^2=1$ is used for `BETADAMP`, as `BETADAMP`$=\dfrac{\alpha_m}{\alpha_d}$.

The pinning files is the same, changing `BETADAMP` to `R0`, where `R0`$=a_0$.

The input file `input/input.in' gives a variety of options. I highly recommend just messing around with to get what each one does.

A output folder is created `out` and inside you will find the average velocity file, and a folder that separetes the average velocity by `BETADAMP` value. You will also find a `positions` folder, that will have the positions files for the particles, the files will be present if and only if `WRITE` on the input file is set to `1`. There will be a `simulator_data.out` file inside, that has informations about the simulation. 

The velocity file has the current value, the velocity in the $x$ direction and the velocity in the $y$ direction, in this order.

One more important thing, the program saves snapshots of the configuration for each current step. This is very useful for cases where the computer shutdown, for example, where I run these simulation there is a serious problem with energy, and this helps a lot to not lose everything that already ran. This can also be used for save space, suppose this: you alread ran an simulation, and has the velocity-force curves, but want some trajectories, how would you do it, if you didn´t save the positions? You could rerun everything, but for many systems this could take a lot of time. With this, you can ran everything, and when you need a trajectory, just set the program to `RECOVERY` (putting a `1` on input file) and run it, with this you will get the desired trajectory. With this you save time and space, just remember to set `WRITE` to `1` too, otherwise it will not write the positions files. By the way, the positions files are very heavy, so there is a `Ncut` on input file, which cuts the positions file by that number. I often use `Ncut: 50`, as the positions file will take only $2\%$ of the original, and the positions will still be very clear and smooth. For animations I recommend using something like `5`, `10`, `20`, it deppends on you.

The integration method can be the Runge-Kutta 4 order, Runge-Kutta 2 order and Euler integration, which one must be defined at compilation time, with `DRK4`, `DRK2` and `DEULER` for the respective method.

For compilation use
```
g++ -std=c++11 -O3 -I ./FileParser -I ./Profiler -I ./headers -Dmethod -o program main.cpp
```

TODO:
- [X] Save/Load system, for backup and use less space.
- [X] Output simulation data.
- [X] Make GSA works.
- [ ] Multithread Simulation.
- [ ] Multithread GSA.
- [ ] Find GSA parameters by files insted of hard coding.
- [ ] Make a system to run many GSA with different parameters for find the absolute ground state.
- [ ] Add references to `README.md`.
- [ ] Add some figures of velocity-force curves, and trajectories.