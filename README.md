# magneto

A fast parallel C++ program for various computations of the 2D [Ising model](http://en.wikipedia.org/wiki/Ising_model). Example uses:

```
magneto -L=30 -TMax=6 -TSteps=10 -en=energy -mag=magnetization
```

For much more examples and things you can do with it, see the **[IPython Notebook](http://nbviewer.ipython.org/github/s9w/magneto/blob/master/physics.ipynb)**.

Needs C++11, OpenMP and boost. Compile with `make`.

## Motivation
There is already a lot of Ising code out there, but most are sub-optimal for various reasons. Among those are slow performance, unreadable code, nonexistant documentation, missing normalization of the physics or missing observables. I hope this program does [slightly better](https://imgs.xkcd.com/comics/standards.png).

## Usage
The **grid size** can be set with `-L=32`, the coupling constant J with `-J=1`, and the number of threads with `-threads=3` (default=3). For performance reasons, `J` is parsed and used as an integer! But usually it's 1 or -1 anyways.

The **temperatures** can be controlled with the three parameters `-TMin=0.1 -TMax=4.53 -TSteps=9`. The odd default `TMax` value corresponds to twice the critical temperature. The interval is divided into `TSteps` steps. Note that start and endpoint are included! It's equivalent to NumPy's `np.linspace(0.0, 3.0, num=4)`.

Besides the default uniform temperature distribution, they can also be generated with a normal distribution around the critical temperature. That helps to concentrate the computation around the interesting point. Use with `-dist=normal`. Minimum and maximum temperatures are always included, the space between is filled with a normal distribution with sigma=1.0.

The **initial state** can be set with `-initial=random`. Valid parameters are `random`, `uniform` (all +1), `neel` (antiferromagnetic or "checkerboard") or a filename. If a filename is given, that configuration is loaded. Important: set the grid length first (`-L=50`)! The format is expected to be ones and zeros with comma seperation. Example:

    1,0,1
    0,1,0
    1,0,1

The algorithm to use in the **equilibration phase** can be set with `-alg1=metro`. Valid parameters are `metro` (Metropolis) and `sw` (Swendsen-Wang). The number of times that algorithm is run during equilibration can be controlled with `-N1=50`.

After the equilibration phase, the algorithm and the number of steps for the **main phase** can be set with `-alg2=metro` and `-N2=500` (works like above). The number of algorithm runs between each step is `-N3=5`.

## Measurements
The parameter `-record=...` controls when and how often measurements are written. By default (`end`), that's once at the end of the main phase. But they can also be recorded as a series during the main phase with `-record=main`. The number of written results is then `N2`.

Every measurement is written into a separate file. On each line, the output files contain the temperature followed by one or more results of the measurement at that temperature. The file is in csv format (=comma separated).

Energy and magnetization can be recorded with `-en=filename` and `-mag=filename` respectively. Those two are properties that can be measured on the system directly and are averaged during the main phase. With `-record=main`, the ith result is the average of all the previous measurements and therefore only the last value is averaged over `N2`. For the default `-record=end` setting, the result is just the average over the `N2` measurements.

**Heat capacity** and **magnetic susceptibility** can be recorded with `-cv=filename` and `-chi=filename` respectively. They are measured from the variation in energy/magnetization (divided by the respective T and T^2 factors).
 
 The **correlation function** can be recorded with `-corr=filename`. It's different in that it outputs not one, but L/2 values. They can be fitted externally to calculate the correlation length and the susceptibility from it. See the notebook for that.

System states can be recorded with `-states=filename`. The output is handled slightly different: A file is written for every temperature like `filename0.txt`, `filename1.txt` and so on. In the `-record=end` mode, the file format is like the initial state format:

    1,0,1
    0,1,0
    1,0,1

With `-record=main`, the different states are written directly under each other. **Example**: Computation is started with `magneto -TMin=2.5 -TSteps=1 -L=80 -measure=states -o=states`. To plot with matplotlib: `plt.imshow(np.loadtxt("states0.txt",  delimiter=","), cmap=plt.cm.Greys, interpolation ="none")`. Result:

![Example state](http://i.imgur.com/xXkFltH.png)