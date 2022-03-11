# Overdamped Particles Simulation With Pinnings Arrays

This code simulates overdamped particles that follows the Langevin equation:

<p align="center"><img src="svgs/5dfc1ce2e62c54d06b56de4743a03846.svg?invert_in_darkmode" align=middle width=320.49067379999997pt height=39.26959575pt/></p>

Where 
<img src="svgs/2953e8a6deff33b641b26972d1dcd6f7.svg?invert_in_darkmode" align=middle width=22.847224949999998pt height=14.611878600000017pt/> 
is the particle velocity, 
<img src="svgs/a2b05ad82c1c76f19c95f9a49d3c8967.svg?invert_in_darkmode" align=middle width=25.5779337pt height=14.15524440000002pt/> 
is the damping constant, that arises from
dissipation processes. 
<img src="svgs/d3aeb8efc0ed9dac5b81dca50ae71a5f.svg?invert_in_darkmode" align=middle width=30.3997056pt height=14.15524440000002pt/> 
is the Magnus term, very important for some overdamped particles, such as magnetic skyrmions.

The right side of the equation are forces: 
<img src="svgs/081c2b57146c2252ef76fb1d4e7c85aa.svg?invert_in_darkmode" align=middle width=40.3915908pt height=27.6567522pt/> 
is the particle-particle interaction,
<img src="svgs/970e831d9bb906eebaacfb804c83ef0e.svg?invert_in_darkmode" align=middle width=37.02932145pt height=31.525041899999984pt/> 
is the particle-pinning interaction and 
<img src="svgs/1732bb70494dab006e4592074fa6d71a.svg?invert_in_darkmode" align=middle width=30.348613349999997pt height=27.6567522pt/> 
is the driven currents force.

For the particle-particle, I modeled the potential as 
<img src="svgs/a8a6daa4de9d144a85b337eeadbb174f.svg?invert_in_darkmode" align=middle width=135.51694695pt height=27.6567522pt/>
, where 
<img src="svgs/e64f3f6406e6eb4787ac528b0405d9dd.svg?invert_in_darkmode" align=middle width=28.732968pt height=22.465723500000017pt/> 
is the modified bessel of first kind, 
<img src="svgs/77ba9956edd663cd452c547f7cd8ca9a.svg?invert_in_darkmode" align=middle width=25.54745655pt height=22.465723500000017pt/> 
is the potential strength of particle 
<img src="svgs/824bb85104c48ef146b8801edfd7fcdd.svg?invert_in_darkmode" align=middle width=15.929626349999998pt height=21.68300969999999pt/>
, and 
<img src="svgs/c0e74a621178592b1b6585beff9d0dd4.svg?invert_in_darkmode" align=middle width=26.390938199999997pt height=14.15524440000002pt/> 
is the distance between particles 
<img src="svgs/eb01ee8fd51d5c8d30f6bc1d644359d6.svg?invert_in_darkmode" align=middle width=13.882435049999998pt height=21.68300969999999pt/> 
and 
<img src="svgs/824bb85104c48ef146b8801edfd7fcdd.svg?invert_in_darkmode" align=middle width=15.929626349999998pt height=21.68300969999999pt/> 
. From this I get the force 
<img src="svgs/1cddac85489fd0afdae8e03ea1f17cd9.svg?invert_in_darkmode" align=middle width=237.50835734999998pt height=27.6567522pt/> 
where 
<img src="svgs/ca9dc377fb43ec08e9a73f657f923c84.svg?invert_in_darkmode" align=middle width=28.732968pt height=22.465723500000017pt/> 
is the modified bessel of second kind.

For the particle-pinning I choose a gaussian potential, with the form 
<img src="svgs/eaeb3397a96b1019b2a0c8a8429e88f7.svg?invert_in_darkmode" align=middle width=159.0478494pt height=57.53473439999999pt/> 
where 
<img src="svgs/c0e74a621178592b1b6585beff9d0dd4.svg?invert_in_darkmode" align=middle width=26.390938199999997pt height=14.15524440000002pt/> 
is the distance between particle 
<img src="svgs/eb01ee8fd51d5c8d30f6bc1d644359d6.svg?invert_in_darkmode" align=middle width=13.882435049999998pt height=21.68300969999999pt/> 
and pinning 
<img src="svgs/824bb85104c48ef146b8801edfd7fcdd.svg?invert_in_darkmode" align=middle width=15.929626349999998pt height=21.68300969999999pt/>
, 
<img src="svgs/7e7b22628f68b36e4f390a8f49f8c84d.svg?invert_in_darkmode" align=middle width=25.995494249999997pt height=22.465723500000017pt/> 
is the pinning strength and 
<img src="svgs/a6fe0988d5149a0280a5db51e58102e2.svg?invert_in_darkmode" align=middle width=23.46090945pt height=14.15524440000002pt/> 
is the pinning radius. From the potential I get the force 
<img src="svgs/b258a08ea8106cdecc724a8fe4cb1d62.svg?invert_in_darkmode" align=middle width=278.7573195pt height=57.53473439999999pt/> 
where 
<img src="svgs/ed798d20c6889223777b27ebde7d667e.svg?invert_in_darkmode" align=middle width=76.87209585pt height=44.70706679999999pt/>.

The last term, is the current force, which can be of any form, in special here I consider the most general form of 
<img src="svgs/26343bc0bcb90a171db6bc0584ea6137.svg?invert_in_darkmode" align=middle width=139.2480738pt height=27.6567522pt/>
. For more detail see the input file `input/input.in`, where all the options are avaiable.

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
lower values less repulsive (if `U0<0` then the particles will be attractive). `BETADAMP` is the parameter from the Magnus force. This is very useful for Skyrmions, for example, if your overdamped particle doesn't have Magnus force, just put `0.0`. The `positions` is the 2D positions for your particles. If you want particles of different species, just put another file with the new species data, and the program will simulate both interacting. One more thing about `BETADAMP`, the normalization 
<img src="svgs/cff46b6fea1eaeb4608893a223469cda.svg?invert_in_darkmode" align=middle width=99.6302934pt height=26.76175259999998pt/> 
is used for `BETADAMP`, as `BETADAMP`
<img src="svgs/4381a894c2903e76be2f2e88a42df875.svg?invert_in_darkmode" align=middle width=42.32654084999999pt height=36.3965877pt/>.

The pinning files is the same, changing `BETADAMP` to `R0`, where `R0`<img src="svgs/17d8b0c8defe55d3706e7b9fdf06aab0.svg?invert_in_darkmode" align=middle width=32.59323209999999pt height=14.15524440000002pt/>.

The input file `input/input.in' gives a variety of options. I highly recommend just messing around with to get what each one does.

A output folder is created `out` and inside you will find the average velocity file, and a folder that separetes the average velocity by `BETADAMP` value. You will also find a `positions` folder, that will have the positions files for the particles, the files will be present if and only if `WRITE` on the input file is set to `1`. There will be a `simulator_data.out` file inside, that has informations about the simulation. 

The velocity file has the current value, the velocity in the <img src="svgs/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode" align=middle width=9.39498779999999pt height=14.15524440000002pt/> direction and the velocity in the <img src="svgs/deceeaf6940a8c7a5a02373728002b0f.svg?invert_in_darkmode" align=middle width=8.649225749999989pt height=14.15524440000002pt/> direction, in this order.

One more important thing, the program saves snapshots of the configuration for each current step. This is very useful for cases where the computer shutdown, for example, where I run these simulation there is a serious problem with energy, and this helps a lot to not lose everything that already ran. This can also be used for save space, suppose this: you alread ran an simulation, and has the velocity-force curves, but want some trajectories, how would you do it, if you didn´t save the positions? You could rerun everything, but for many systems this could take a lot of time. With this, you can ran everything, and when you need a trajectory, just set the program to `RECOVERY` (putting a `1` on input file) and run it, with this you will get the desired trajectory. With this you save time and space, just remember to set `WRITE` to `1` too, otherwise it will not write the positions files. By the way, the positions files are very heavy, so there is a `Ncut` on input file, which cuts the positions file by that number. I often use `Ncut: 50`, as the positions file will take only <img src="svgs/45a0b00b513fa74f40b37aafadb94773.svg?invert_in_darkmode" align=middle width=21.91788224999999pt height=24.65753399999998pt/> of the original, and the positions will still be very clear and smooth. For animations I recommend using something like `5`, `10`, `20`, it deppends on you.

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