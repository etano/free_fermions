# free_fermions

Calculate various properties of free fermions

# Building

    apt-get update
    apt-get install build-essential \
                    libmpfrc++-dev \
                    libmpfr-dev \
                    libgsl-dev \
                    libblas-dev \
                    liblapack-dev \
    git clone https://github.com/etano/free_fermions.git
    cd free_fermions/
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

# Output

- sign.dat: average value of the sign

    root@simpimc:~/free_fermions# head sign.dat
    0.93190268

- PpF(B).dat: fermion(boson) permutation sector weights

    root@simpimc:~/free_fermions# head PpF.dat
    0 1.0359151
    1 -0.036530109
    2 0.0003147998
    3 0.00030671038
    4 -2.913327e-06
    5 -3.1717081e-06
    6 2.3761677e-08
    7 -5.1503409e-07
    8 1.4676348e-08
    9 5.4664663e-09

- Pk.dat: cycle length probabilities

    root@simpimc:~/free_fermions# head Pk.dat
    0 0.99498654
    1 0.0049705185
    2 4.2547098e-05
    3 3.9177657e-07
    4 3.1847026e-09
    5 1.9145149e-11
    6 6.1849149e-14

- Pl.dat: permutation length weights (see paper)

    root@simpimc:~/free_fermions# head Pl.dat
    0 210.55151
    1 74.443064
    2 40.52167
    3 26.319597
    4 18.832771
    5 14.326575
    6 11.369006

- EF(B).dat: fermion(boson) total energy

    root@simpimc:~/free_fermions# head EF.dat
    3.8555463

- EF(B)N.dat: fermion(boson) energy per particle

    root@simpimc:~/free_fermions# cat EFN.dat
    0.55079232

- EF(B)_thermo.dat: fermion(boson) thermodynamic limit energy per particle

    root@simpimc:~/free_fermions# cat EF_thermo.dat
    1.1026884

- EpF(B).dat: fermion(boson) permutation sector total energy

    root@simpimc:~/free_fermions# cat EpF.dat
    0 -3.8362006
    1 -3.2882681
    2 -2.7402395
    3 -2.7403358
    4 -2.1922109
    5 -2.1923071
    6 -1.6441822
    7 -2.1924033
    8 -1.6442784
    9 -1.6442784
    10 -1.0961532
    11 -1.6443747
    12 -1.0962498
    13 -1.0962498
    14 -0.5481205

- El.dat: cycle length energy weight (see paper)

    root@simpimc:~/free_fermions# cat El.dat
    0 0.54802866
    1 0.27406244
    2 0.1827083
    3 0.13703122
    4 0.10962498
    5 0.091354092
    6 0.078302928

- PpEpF(B).dat: fermion(boson) permutation sector probability and total energy (same as above, only combined)

    root@simpimc:~/free_fermions# cat PpEpF.dat
    0 1.0359151 -3.8362006
    1 -0.036530109 -3.2882681
    2 0.0003147998 -2.7402395
    3 0.00030671038 -2.7403358
    4 -2.913327e-06 -2.1922109
    5 -3.1717081e-06 -2.1923071
    6 2.3761677e-08 -1.6441822
    7 -5.1503409e-07 -2.1924033
    8 1.4676348e-08 -1.6442784
    9 5.4664663e-09 -1.6442784
    10 -1.430855e-10 -1.0961532
    11 2.6629972e-09 -1.6443747
    12 -3.9901075e-11 -1.0962498
    13 -2.5294813e-11 -1.0962498
    14 4.6224327e-13 -0.5481205
