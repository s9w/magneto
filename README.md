# ising2d

A C++ program for calculating various computations of the 2D ising model. Example uses:

```
make
ising -L=60 -TSteps=9 -avgN=500 -measure=energy -o=ising_e
ising -L=60 -TSteps=12 -avgN=500 -measure=mag -o=ising_m
```

It randomly (1 or -1) initializes the grid with the specified size and each of the temperatures. The equilibration is done with the Swendsen-Wang algorithm.

The options and their defaults are:
- `-L=30`: The size of the grid. Produces grid of size L*L
- `-TMin=0.1`: Minimal temperature
- `-TMin=4.53`: Maximal temperature. Approx 2*T_critical
- `-TSteps=10`: Number of temperature steps. Divides the range into 10 pieces
