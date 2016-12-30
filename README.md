# wegner_ising

This is some Monte Carlo code that can be used on the 3D Wegner Ising model. There are a few measurement kernels available (some a WIP) for use in measuring things like heat capacity or a Wilson loop. The code is modularized to make it easy to add new kernels (uses function pointers that get passed to the MC evolution function).
