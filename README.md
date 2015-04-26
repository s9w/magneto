# ising2d

A fast parallel C++ program for various computations of the 2D [Ising model](http://en.wikipedia.org/wiki/Ising_model). Example uses:

```
ising -L=30 -TSteps=9 -measure=energy -o=energy
ising -L=30 -TSteps=12 -measure=mag -o=magnetization
```

OpenMP is used for the parallel work as well as some parts of C++11. No other external libraries needed.

## Motivation
There is already a lot of Ising code out there, but most are sub-optimal for various reasons. Among those are slow performance, unreadable code, nonexistant documentation, missing normalization of the physics or missing observables. I hope this program does slightly better.

## Usage
Output filename can be set with `-o=abc`. That'll be extended into `abc.txt`.

The grid size can be set with `-L=32`.

The temperatures are controled with three parameters: Min and max temperature and the number of steps. So `-TMin=0.0 -TMax=3.0 -TSteps=4` would do calculate for T=0.0, T=1.0, T=2.0, T=3.0. It's equivalent to numpy's `np.linspace(0.0, 3.0, num=4, endpoint=True)`.

The grid is randomly initialized with 1 and -1 values. The equilibration is done with the Swendsen-Wang algorithm. From there, a number of different computations can be done. Also see the **[IPython Notebook](http://nbviewer.ipython.org/github/s9w/ising2d/blob/master/usage.ipynb)** for examples. The number of threads can be controlled with `-threads=4`.

### States
To output system states, use `-measure=states`. A number of csv files equivalent to `TSteps` will be written. They contain 0 and 1 corresponding to the spin states.

**Example**: Computation is started with `ising -TMin=2.5 -TSteps=1 -L=80 -measure=states -o=states`. To plot with matplotlib: `plt.imshow(np.loadtxt("states0.txt",  delimiter=","), cmap=plt.cm.Greys, interpolation ="none")`. Result:

![Example state](http://i.imgur.com/xXkFltH.png)

### Energy
`-measure=energy` measures the energy of the system given by the hamiltonian. It's averaged over `avgN` measurements of different systems.

### Magnetization
`-measure=mag` measures the absolute value of the magnetization. It's averaged over `avgN` measurements of different systems.

### Heat capacity
`-measure=cv` measures the specific heat capacity from energy variations. Calculated from `avgN` different systems.

### Magnetic Susceptibility
`-measure=cv` measures the specific heat capacity from the variation of the magnetization. Calculated from `avgN` different systems.

## Options
Overview of all parameters and their defaults:
- `-L=30`: The size of the grid. Produces grid of size L*L
- `-TMin=0.1`: Minimal temperature
- `-TMin=4.53`: Maximal temperature. Approx 2*T_critical
- `-TSteps=10`: Number of temperature steps. Divides the range into 10 pieces
- `-avgN=500`: Number of systems to use for computations
- `-o=out`: Filename base for the results. `.txt` gets appended 
- `-threads=2`: Number of parallel threads
- `-measure=energy`: Type of meausurement. Valid parameters are: `energy`, `mag` (magnetization), `cv` (heat capacity), `chi` (magnetic susceptibility), `corrLen` (correlation Length), `states` (system states)