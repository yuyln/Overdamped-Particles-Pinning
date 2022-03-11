# Overdamped Particles Simulation With Pinnings Arrays

This code simulates overdamped particles that follows the Langevin equation:

<p align="center"><img src="svgs/885e61cb3afce16640aae6a0e448d303.svg?invert_in_darkmode" align=middle width=320.49067379999997pt height=39.26959575pt/></p>

Where 
<img src="svgs/5474c8baa2bc6feabb0eac4237772aab.svg?invert_in_darkmode" align=middle width=14.628015599999989pt height=14.611878600000017pt/> 
is the particle velocity, 
<img src="svgs/923122acc2287bfda324405831221148.svg?invert_in_darkmode" align=middle width=17.35872434999999pt height=14.15524440000002pt/> 
is the damping constant, that arises from
dissipation processes. 
<img src="svgs/812375f957f3b78a317a4d7150c3ae73.svg?invert_in_darkmode" align=middle width=22.18049624999999pt height=14.15524440000002pt/> 
is the Magnus term, very important for some overdamped particles, such as magnetic skyrmions.

The right side of the equation are forces: 
<img src="svgs/e95e66507f099399896b77ab8f29283d.svg?invert_in_darkmode" align=middle width=32.17238144999999pt height=27.6567522pt/> 
is the particle-particle interaction,
<img src="svgs/782afeb48b808fd9ab3f9c49d375d7f8.svg?invert_in_darkmode" align=middle width=28.810112099999987pt height=31.525041899999984pt/> 
is the particle-pinning interaction and 
<img src="svgs/c607da508083ad5701363425633a95f0.svg?invert_in_darkmode" align=middle width=22.12940399999999pt height=27.6567522pt/> 
is the driven currents force.

For the particle-particle, I modeled the potential as 
<img src="svgs/e3e70515565a7baf69520fdd8ea67d9e.svg?invert_in_darkmode" align=middle width=127.29773759999999pt height=27.6567522pt/>
, where 
<img src="svgs/e9aaf76f0a4e315c1f4a8f8cab730a2d.svg?invert_in_darkmode" align=middle width=20.513758649999993pt height=22.465723500000017pt/> 
is the modified bessel of first kind, 
<img src="svgs/22e03b2e0480e5d934a24177ee642c6e.svg?invert_in_darkmode" align=middle width=17.32824719999999pt height=22.465723500000017pt/> 
is the potential strength of particle 
<img src="svgs/36b5afebdba34564d884d347484ac0c7.svg?invert_in_darkmode" align=middle width=7.710416999999989pt height=21.68300969999999pt/>
, and 
<img src="svgs/92e0822b1528090efc2435d2ae60c9ee.svg?invert_in_darkmode" align=middle width=18.17172884999999pt height=14.15524440000002pt/> 
is the distance between particles 
<img src="svgs/77a3b857d53fb44e33b53e4c8b68351a.svg?invert_in_darkmode" align=middle width=5.663225699999989pt height=21.68300969999999pt/> 
and 
<img src="svgs/36b5afebdba34564d884d347484ac0c7.svg?invert_in_darkmode" align=middle width=7.710416999999989pt height=21.68300969999999pt/> 
. From this I get the force 
<img src="svgs/3ee4eeebf0289c1b1309ab62599acf79.svg?invert_in_darkmode" align=middle width=229.28914799999998pt height=27.6567522pt/> 
where 
<img src="svgs/6d6968ce16fa7a6d8af6261f0d09b3b9.svg?invert_in_darkmode" align=middle width=20.513758649999993pt height=22.465723500000017pt/> 
is the modified bessel of second kind.

For the particle-pinning I choose a gaussian potential, with the form 
<img src="svgs/660a9f214c81aaa128108051f54a9868.svg?invert_in_darkmode" align=middle width=150.82864005pt height=57.53473439999999pt/> 
where 
<img src="svgs/92e0822b1528090efc2435d2ae60c9ee.svg?invert_in_darkmode" align=middle width=18.17172884999999pt height=14.15524440000002pt/> 
is the distance between particle 
<img src="svgs/77a3b857d53fb44e33b53e4c8b68351a.svg?invert_in_darkmode" align=middle width=5.663225699999989pt height=21.68300969999999pt/> 
and pinning 
<img src="svgs/36b5afebdba34564d884d347484ac0c7.svg?invert_in_darkmode" align=middle width=7.710416999999989pt height=21.68300969999999pt/>
, 
<img src="svgs/1049ded9e1be03670fc3963966339893.svg?invert_in_darkmode" align=middle width=17.77628489999999pt height=22.465723500000017pt/> 
is the pinning strength and 
<img src="svgs/007094eee0f16d09ce121fc2ba8e7107.svg?invert_in_darkmode" align=middle width=15.24170009999999pt height=14.15524440000002pt/> 
is the pinning radius. From the potential I get the force 
<img src="svgs/9a8051ff565cb120a39e222aa144fe54.svg?invert_in_darkmode" align=middle width=270.53811014999997pt height=57.53473439999999pt/> 
where 
<img src="svgs/66ea7bb93062b841d8ee79b99856202e.svg?invert_in_darkmode" align=middle width=68.65288649999998pt height=44.70706679999999pt/>.

The last term, is the current force, which can be of any form, in special here I consider the most general form of 
<img src="svgs/a8181d56bbd2cd2afe6e379912da1149.svg?invert_in_darkmode" align=middle width=131.02886445pt height=27.6567522pt/>
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
<img src="svgs/8b1a0facd26bc55bbe29cb938fb04c46.svg?invert_in_darkmode" align=middle width=91.41108404999999pt height=26.76175259999998pt/> 
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