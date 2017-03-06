# GonioTests IDL Library

## Disclaimer

This code has been tested with [IDL (Interactive Data Language)](http://www.harrisgeospatial.com/ProductsandSolutions/GeospatialProducts/IDL.aspx), 
a commercial data processing language. GDL (GNU Data Language) is an opensource version of IDL. 
This code runs with GDL (tested on version 0.9.6).

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
 
### Setting a set of antenna for your simulation

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
% Compiled module: ANTENNA_SET__DEFINE.
% Compiled module: ANTENNA__DEFINE.
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
IDL> write_antenna_set,file='test_antenna',ant_set=ant,path='data/temp'
% Compiled module: WRITE_ANTENNA_SET.
% WRITE_ANTENNA_SET: You entered the followind parameters : 
% WRITE_ANTENNA_SET: X+ antenna : h=1.00 al= 90.0 be= 45.0
% WRITE_ANTENNA_SET: X- antenna : h=1.00 al= 90.0 be=-45.0
% WRITE_ANTENNA_SET: Z  antenna : h=1.00 al=  0.0 be=  0.0
% WRITE_ANTENNA_SET: Dipole     : h=1.00 al= 90.0 be= 90.0
% WRITE_ANTENNA_SET: Antenna set parameters writen in data/temp/test_antenna.ant
```
This sequence of instructions will create an antenna file with the desired 
parameters at `data/temp/test_antennna.ant`.

Note that when launched with no parameters, the user can enter the various 
parameters manually.

### Simulating goniopolarimetric data 

The main routine is `dfb_test_run`. A minimal set of keywords to specify is 
shown in the example below. The source code of the script gives all details 
on the various keywords and options.

```idl
IDL> dfb_test_run,antenna_file='test_antenna',file_ext='dummy',output_path='data/'
% Compiled module: DFB_TEST_RUN.
% Compiled module: QUATERNION.
The quaternion.pro procedure set has been successfully compiled and loaded ! 
 (c) BC, Jun 11 2002
Polar :          208 -- Angles :         2522 -- Flux :            2
 -> Total :      1049152
% Compiled module: READ_ANTENNA_SET.
% Compiled module: ANTENNA_SET__DEFINE.
% Compiled module: ANTENNA__DEFINE.
% READ_ANTENNA_SET: Reading Antenna Parameters from : data/temp//test_antenna.ant
% READ_ANTENNA_SET:  Converting angles into radian ...
% READ_ANTENNA_SET: Antenna Parameters loaded.
% Compiled module: DATA_N2__DEFINE.
% Compiled module: DATA_N3B__DEFINE.
% Compiled module: DATA_EPHEM__DEFINE.
% Compiled module: V1V2_XTND.
% Compiled module: V1IV2_XTND.
% Compiled module: WRITE_DATA_BINARY.
% Compiled module: WRITE_ANTENNA_SET.
% WRITE_ANTENNA_SET: You entered the followind parameters : 
% WRITE_ANTENNA_SET: X+ antenna : h=1.00 al= 90.0 be= 45.0
% WRITE_ANTENNA_SET: X- antenna : h=1.00 al= 90.0 be=-45.0
% WRITE_ANTENNA_SET: Z  antenna : h=1.00 al=  0.0 be=  0.0
% WRITE_ANTENNA_SET: Dipole     : h=1.00 al= 90.0 be= 90.0
% WRITE_ANTENNA_SET: Antenna set parameters writen in data/temp/dfb_test_dummy.ant
% Compiled module: DFB_MAIN.
% READ_ANTENNA_SET: Reading Antenna Parameters from : data/temp//test_antenna.ant
% READ_ANTENNA_SET:  Converting angles into radian ...
% READ_ANTENNA_SET: Antenna Parameters loaded.
% Compiled module: MAKE_VECT_SPH.
% Compiled module: CROSSP1.
% Compiled module: ANGULAR_DISTANCE.
% Program caused arithmetic error: Floating illegal operand
IDL> 
```
The script used the antenna file prepared in the previous section, and wrote 
the following files:

* `data/n2/Pdfb_test_dummy.00` containing the *level 2* modeled measurements
* `data/n3b/N3b_Ixx_test_dummy.00` containing the input wave characteristics
* `data/n3b/N3b_Oxx_test_dummy.00` containing the reconstructed wave characteristics
* `data/ephem/dfb_test_dummy.ephem` containing the input radio source location

### Plotting results

Currently, there is no integrated routine to plot data, but an example script
is provided to help you to plot out the resulting data. The script 
`juice_rwi_boom.pro` is running under IDL and in GDL (after a dirty
fix in the `read_data_binary.pro` routine). This script shows how to run 
simulation, load data and plot data. 
