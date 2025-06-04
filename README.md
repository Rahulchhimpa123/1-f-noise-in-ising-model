This code was used in the Physical Review E paper titled "$1/f$ noise in the Ising model" by Rahul Chhimpa and Avinash Chand Yadav.
If you use this code in your work, please cite the above article. 
All of the simulations are performed in FORTRAN, while the visualization is in Python. We provide Fortran codes, Python files, and data files.
------------------------------------------------------------------------------------------------------------------------\

In Fortran codes, the intrinsic routine generates random numbers using the `xoshiro256 pseudorandom number generator (PRNG)`. 

-------------------------------------------------------------------------------------------------------------------------
Folder -- Initial_configuration

This Folder contains a Fortran code `initial_configuration.f90`, which will generate two files named `fort.10` and `fort.11`. By running the Python code `fig_ising_conf_1.py`, we can read the data file and generate the `Figure 1`, used in the above article.

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
Folder -- Pdf

This folder contains the `distribution.f90` code, which has a subroutine to compute the probability distribution of the data. This subroutine can be easily paired with the `initial_configuration.f90` file to read the array containing a number of spin flips per Monte Carlo step data and produce probability distribution. It is used to generate `Figure_2` in the article.

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
Folder -- correlation

This folder contains the `correlation.f90` code, which is used to find the two-time autocorrelation function of the data array. Ensemble averaging is also performed and parallel programming is used by using `OpenMP`. We also provide data files and Python code to plot the correlation function and perform necessary scaling operations.


-----------------------------------------------------------------------------------------------------------------------------------
Folder -- psd

This folder contains a subroutine to find the power spectral density (PSD) of the signal. We also provide the Python file to plot PSD and perform the scaling to get the data collapse.
