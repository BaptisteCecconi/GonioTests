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

### Setting an set of antenna for your simulation

The antenna files are binary files containing 
a series of four records composed of three 32 bit floating point values (little 
endian) corresponding to the length, latitude and azimuth of the antenna 
vectors. Only the first three antenna vectors of the antenna file are used with 
the provided set of routines. The antenna are named `Xp`, `Xm`, `Z` and `Dip`, 
which are the names of the antenna of the Cassini/RPWS instrument. 

The `write_antenna_set` routine can be used to produce antenna (`.ant`) files 
containing antenna parameters. When launched with no parameters, the user can 
enter the various parameters manually:

```idl
IDL> write_antenna_set
```

### Simulating goniopolarimetric data 

The main routine is `dfb_test_run`.