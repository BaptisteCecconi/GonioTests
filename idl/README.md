# GonioTests IDL Library

## Disclaimer

This code has been tested with [IDL (Interactive Data Language)](http://www.harrisgeospatial.com/ProductsandSolutions/GeospatialProducts/IDL.aspx), 
a commercial data processing language. GDL (GNU Data Language) is an opensource version of IDL. This code should work with GDL.

## Usage

This library can be used to evaluate the reconstruction of incoming waves 
characteristics (flux, polarization, direction of arrival, planarity), when
measured with a goniopolarimetric radio instrument. The input radio wave 
characteristics are configured as a series of Stokes parameters and direction
of propagation. The instrument is defined by the effective sensor directions 
and lengths. The following noises and biases are configurable:
* Error on each antenna length calibration
* Error on each antenna direction calibration
* Digitization effect for 8 and 12 floating point compression
* Measurement noise with selectable spectral bandwidth and integration time
* Change in the wave flux between successive measurements
* Error in the background level determination
* Error on the cross-channel phase calibration
* Non perfect dipole null

The frequency of observation is not a parameter, as we consider observations in
the short dipole range where the antenna response is constant over frequencies.

### Data directories

The library is working smoothly if a specific data directory structure is 
followed. A `data/` directory must be created. This will be the main directory
for input and output files. It must contain four subdirectories:

* `data/n2/`  contains the simulated goniopolarimetric data, which are auto- 
and cross-correlations (*level 2 data*).  
* `data/n3b/` contains the simulated reconstructed data, which are Stokes 
parameters and direction of arrival using `DFb`inversion (*level 3b data*).  
* `data/ephem/` contains the input waves direction of propagations.  
* `data/temp/` contains the test antenna parameter files. 

The data are stored in a custom fixed-length record binary format. 
 
### Setting an set of antenna for your simulation

The antenna files are binary files containing 
a series of four records composed of three 32 bit floating point values (little 
endian) corresponding to the length, latitude and azimuth of the antenna 
vectors. Only the first three antenna vectors of the antenna file are used with 
the provided set of routines. The antenna are named `Xp`, `Xm`, `Z` and `Dip`, 
which are the names of the antenna of the Cassini/RPWS instrument. 

The `write_antenna_set` routine can be used to produce antenna (`.ant`) files 
containing antenna parameters. 

```idl
IDL> ant = {antenna_set}
IDL> ant.xp.h = 1.0
IDL> ant.xp.al = 90.0
IDL> ant.xp.be = 45.0
IDL> ant.xm.h = 1.0
IDL> ant.xm.al = 90.0
IDL> ant.xm.be = -45.0
IDL> ant.z.h = 1.0
IDL> ant.z.al = 0.0
IDL> ant.z.be = 0.0
IDL> ant.dip.h = 1.0
IDL> ant.dip.al = 90.0
IDL> ant.dip.be = 90.0
IDL> write_antenna_set,file='test_antenna',ant_set=ant,path='data/temp/'
```
This sequence of instructions will create an antenna file with the desired 
parameters at `data/temp/test_antennna.ant`.

Note that when launched with no parameters, the user can enter the various 
parameters manually.

### Simulating goniopolarimetric data 

The main routine is `dfb_test_run`.