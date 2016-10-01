# free_fermions
Calculate various properties of free fermions

# Building

    make

# Parameters

All command line parameters have the form:

    -PARAMETER VALUE

Available parameters are:

    -N : INT : Number of particles
    -D : INT : Physical dimension
    -units : BOOL : Optional. Whether or not to use electron gas units.
    -n_max : INT : Optional. Number of periodic images to include, default 10.
    -precision : INT : Optional. Number of decimal places to include, default 8.
    -rs : DOUBLE : If "-units" is set, then this is the Wigner-Seitz radius
    -theta : DOUBLE : If "-units" is set, then this is T/T_F.
    -L : DOUBLE : If "-units" is not set, then this is the box side length.
    -T : DOUBLE : If "-units" is not set, then this is the temperature.
    -lambda : DOUBLE : If "-units" is not set, then this is hbar^2/2m.

Example:

    ./free_fermions -N 7 -D 3 -polarized 1  -L 10.0 -units 1 -rs 8.0 -theta 8.0
