# Overdamped Particles Simulation With Pinnings Arrays

This code simulates overdamped particles that follows the Langevin equation:

<p align="center"><img src="svgs/8247f6caa7b04141eed441d42a66cbbf.svg?invert_in_darkmode" align=middle width=362.8324095pt height=39.26959575pt/></p>

Where 
<img src="svgs/0374c41447147ca3648ee247708ef2dc.svg?invert_in_darkmode" align=middle width=22.847224949999998pt height=14.611878600000017pt/> 
is the particle velocity, 
<img src="svgs/2c81a6cce979f48f11d65a4ff4c3e36a.svg?invert_in_darkmode" align=middle width=25.5779337pt height=14.15524440000002pt/> 
is the damping constant, that arises from
dissipation processes. 
<img src="svgs/b76b989c74373553c3c50388d3eebfe0.svg?invert_in_darkmode" align=middle width=30.3997056pt height=14.15524440000002pt/> 
is the Magnus term, very important for some overdamped particles, such as magnetic skyrmions.

The right side of the equation are forces: 
<img src="svgs/dc3bdd931bca1c2d6c2fc539a722a541.svg?invert_in_darkmode" align=middle width=40.3915908pt height=27.6567522pt/> 
is the particle-particle interaction,
<img src="svgs/f904f960c240cf947b32f1d0c6c34e24.svg?invert_in_darkmode" align=middle width=37.02932145pt height=31.525041899999984pt/> 
is the particle-pinning interaction and 
<img src="svgs/e5234634c204183fce71ddd095d239a4.svg?invert_in_darkmode" align=middle width=30.348613349999997pt height=27.6567522pt/> 
is the driven currents force.

For the particle-particle, I modeled the potential as 
<img src="svgs/e316b940aa43ecae7e77ffbb3042f466.svg?invert_in_darkmode" align=middle width=135.51694695pt height=27.6567522pt/>
, where 
<img src="svgs/24c5603126b66cc1f9a1f7dc68f456b7.svg?invert_in_darkmode" align=middle width=28.732968pt height=22.465723500000017pt/> 
is the modified bessel of first kind, 
<img src="svgs/c1ad03439063568075616325b8262b1c.svg?invert_in_darkmode" align=middle width=25.54745655pt height=22.465723500000017pt/> 
is the potential strength of particle 
<img src="svgs/317a2c08d9f6e7dd873fc65fbee394e3.svg?invert_in_darkmode" align=middle width=15.929626349999998pt height=21.68300969999999pt/>
, and 
<img src="svgs/fcdb70efae8e21c43e174377b6f4c03d.svg?invert_in_darkmode" align=middle width=26.390938199999997pt height=14.15524440000002pt/> 
is the distance between particles 
<img src="svgs/14afed0aec1ac4b185b819a4b510bff6.svg?invert_in_darkmode" align=middle width=13.882435049999998pt height=21.68300969999999pt/> 
and 
<img src="svgs/317a2c08d9f6e7dd873fc65fbee394e3.svg?invert_in_darkmode" align=middle width=15.929626349999998pt height=21.68300969999999pt/> 
. From this I get the force 
<img src="svgs/baf9d687c9ec72993b3a331d0fe14a9b.svg?invert_in_darkmode" align=middle width=237.50835734999998pt height=27.6567522pt/> 
where 
<img src="svgs/ad2051ae86a51ad37b859e099fe609d7.svg?invert_in_darkmode" align=middle width=28.732968pt height=22.465723500000017pt/> 
is the modified bessel of second kind.

For the particle-pinning I choose a gaussian potential, with the form 
<img src="svgs/9cbf620c7582b40367ca1c4cd44d1378.svg?invert_in_darkmode" align=middle width=171.8332836pt height=57.53473439999999pt/> 
where 
<img src="svgs/fcdb70efae8e21c43e174377b6f4c03d.svg?invert_in_darkmode" align=middle width=26.390938199999997pt height=14.15524440000002pt/> 
is the distance between particle 
<img src="svgs/14afed0aec1ac4b185b819a4b510bff6.svg?invert_in_darkmode" align=middle width=13.882435049999998pt height=21.68300969999999pt/> 
and pinning 
<img src="svgs/317a2c08d9f6e7dd873fc65fbee394e3.svg?invert_in_darkmode" align=middle width=15.929626349999998pt height=21.68300969999999pt/>
, 
<img src="svgs/9387b5548255b313ba3bdb0e512b3e1c.svg?invert_in_darkmode" align=middle width=25.995494249999997pt height=22.465723500000017pt/> 
is the pinning strength and 
<img src="svgs/93453c892be73e29b3db9df48b14587b.svg?invert_in_darkmode" align=middle width=23.46090945pt height=14.15524440000002pt/> 
is the pinning radius. From the potential I get the force 
<img src="svgs/ce2257caadd72882f246be515eaa35ab.svg?invert_in_darkmode" align=middle width=291.5427471pt height=57.53473439999999pt/> 
where 
<img src="svgs/e7e8964f640e4fbcbe0423b8a1a424cc.svg?invert_in_darkmode" align=middle width=76.87209585pt height=44.70706679999999pt/>.

The third term, is the current force, which can be of any form, in special here I consider the most general form of 
<img src="svgs/cd53890f1d3199b999128b63de2fc5de.svg?invert_in_darkmode" align=middle width=139.2480738pt height=27.6567522pt/>
. For more detail see the input file `input/input.in`, where all the options are avaiable.

The last term is the temperature effect, it is modeled using random kicks in the particle position, and we use a distribuition such that <img src="svgs/41d1475bd1f4d4c7c638e2034f307782.svg?invert_in_darkmode" align=middle width=77.3644806pt height=27.94539330000001pt/>.

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
<img src="svgs/d289beee1f3b63d5c3257e928343fb3f.svg?invert_in_darkmode" align=middle width=99.6302934pt height=26.76175259999998pt/> 
is used for `BETADAMP`, as `BETADAMP`
<img src="svgs/c6412ef4a1bf3c4ea666b2bee45af636.svg?invert_in_darkmode" align=middle width=50.54575019999999pt height=36.3965877pt/>.

The pinning files is the same, changing `BETADAMP` to `R0`, where `R0`<img src="svgs/3fc6810807d36c390b212098f3e5222d.svg?invert_in_darkmode" align=middle width=40.812441449999994pt height=14.15524440000002pt/>.

The input file `input/input.in' gives a variety of options. I highly recommend just messing around with to get what each one does.

A output folder is created `out` and inside you will find the average velocity file, and a folder that separetes the average velocity by `BETADAMP` value. You will also find a `positions` folder, that will have the positions files for the particles, the files will be present if and only if `WRITE` on the input file is set to `1`. There will be a `simulator_data.out` file inside, that has informations about the simulation. 

The velocity file has the current value, the velocity in the <img src="svgs/cd0781ba124922cc88b8210a09b4162e.svg?invert_in_darkmode" align=middle width=17.61419715pt height=14.15524440000002pt/> direction and the velocity in the <img src="svgs/7d226986593183423b49f90d8303043c.svg?invert_in_darkmode" align=middle width=16.8684351pt height=14.15524440000002pt/> direction, in this order.

One more important thing, the program saves snapshots of the configuration for each current step. This is very useful for cases where the computer shutdown, for example, where I run these simulation there is a serious problem with energy, and this helps a lot to not lose everything that already ran. This can also be used for save space, suppose this: you alread ran an simulation, and has the velocity-force curves, but want some trajectories, how would you do it, if you didn´t save the positions? You could rerun everything, but for many systems this could take a lot of time. With this, you can ran everything, and when you need a trajectory, just set the program to `RECOVERY` (putting a `1` on input file) and run it, with this you will get the desired trajectory. With this you save time and space, just remember to set `WRITE` to `1` too, otherwise it will not write the positions files. By the way, the positions files are very heavy, so there is a `Ncut` on input file, which cuts the positions file by that number. I often use `Ncut: 50`, as the positions file will take only <img src="svgs/45a0b00b513fa74f40b37aafadb94773.svg?invert_in_darkmode" align=middle width=21.91788224999999pt height=24.65753399999998pt/> of the original, and the positions will still be very clear and smooth. For animations I recommend using something like `5`, `10`, `20`, it deppends on you.

The integration method can be the Runge-Kutta 4 order, Runge-Kutta 2 order and Euler integration, which one must be defined at compilation time, with `-DRK4`, `-DRK2` and `-DEULER` for the respective method.

One last thing, at any time you can modify the functions of interactions, you can do that by changing the lambdas called when creating the `Simulator` object. Pay attention that there is only 3 boxes, one for pinning, one for particle forces and one for particle potential. This is because the pinning potential and force are very similar, and therefore needing only one box scheme. If this is not your case, you will need to change the code by adding one more box and changing on `Force` and `Potential` functions.

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
- [ ] Find GSA parameters by files instead of hard coding.
- [ ] Make a system to run many GSA with different parameters for find the absolute ground state.
- [ ] Add references to `README.md`.
- [ ] Add some figures of velocity-force curves, and trajectories.