biHelix
=======
A program implemented in Fortran 95/2003 is to find the properties of an alpha 
helix of a protein.

The algorithm proposed by Sugeta et al is employed. This implementation 
refers to helanal (in Fortran 77) developed by Kumar et al. 
Except for the programming language, the main difference to helanal is that 
Singular Value Decomposition is employed other than matrix inverse in helanal
to find the helical axis.

Reference:
 Kumar, S. and Bansal, M. (1996). Structural and Sequence Characteristics
of Long Alpha Helices in Globular Proteins. Biophysical J.,71, 1574-1586.

 Sugeta H, Miyazawa T. General method for calculating helical parameters of
polymer chains from bond lengths, bond angles, and internal-rotation angles.
Biopolymers. 1967;5: 673–679.

Prerequisites
-------------

1. cmake
2. Lapack
3. A fortran compiler (gfortran v5.3 had been tested.)

Compilation
------------
    cmake ../biHex
    make

The executable file will be generated and saved to the bin directory;

Command line arguments
-----------------------
General usage:

    biHelix [-i=required] [-o=helix_out.txt] [-q] [-v]

 Options and flags {default values}:  
 -i                     --in:    The input file name  {required}  
 -o                    --out:    The output file name  {helix_out.txt}  
 -q                  --quiet:    Only output to file  
 -v                --verbose:    write more to file  

 Positional arguments:  
 None  

 Also, -?, -h, -H, -help, --help, and --usage are recognized.  

 Input format
 -------------
The first line of an input file should contain two numbers, the number of atoms and "3". The following lines obey the ATOM lines of PDB format.
The following example as an input file is the segment which contains CA atoms of residue 24-36 of chain C from 1dxr.pdb.  

13   3  
ATOM    171 CA   HIS C  24      90.466  50.611   7.342  
ATOM    181 CA   PRO C  25      86.705  50.551   7.382  
ATOM    188 CA   ALA C  26      86.402  48.355   4.259  
ATOM    193 CA   THR C  27      88.767  45.715   5.625  
ATOM    200 CA   VAL C  28      86.811  45.515   8.869  
ATOM    207 CA   LYS C  29      83.557  45.184   6.898  
ATOM    216 CA   ALA C  30      84.997  42.387   4.777  
ATOM    221 CA   LYS C  31      86.121  40.522   7.889  
ATOM    230 CA   LYS C  32      82.787  41.128   9.641  
ATOM    239 CA   GLU C  33      81.053  39.523   6.681  
ATOM    248 CA   ARG C  34      83.361  36.501   6.568  
ATOM    259 CA   ASP C  35      82.888  36.000  10.284  
ATOM    267 CA   ALA C  36      79.146  36.539  10.043  

Output format
--------------
Directions:
  1       -0.412        -0.898         0.156  
  2       -0.429        -0.888         0.167  
  3       -0.443        -0.887         0.133  
  4       -0.457        -0.879         0.136  
  5       -0.458        -0.880         0.129  
  6       -0.509        -0.852         0.121  
  7       -0.522        -0.826         0.210  
  8       -0.468        -0.851         0.240  
  9       -0.462        -0.857         0.229  
 10       -0.501        -0.844         0.192  
 Helix origins:  
  1       88.223        49.613         5.994  
  2       87.569        48.161         6.222  
  3       86.889        46.776         6.441  
  4       86.262        45.460         6.664  
  5       85.581        44.121         6.833  
  6       84.857        42.801         7.100  
  7       84.258        41.555         7.321  
  8       83.471        40.182         7.626  
  9       82.860        38.889         7.954  
 10       82.137        37.611         8.254  
 11       81.415        36.469         8.500  
 Reference axis:  
       0.000         0.000         1.000  
Tilt angle w.r.t Reference axis:      80.385  
ATOM    171 CA   HIS C  24      90.466  50.611   7.342   0.000   0.000   0.000 
ATOM    181 CA   PRO C  25      86.705  50.551   7.382   0.000   0.000   0.000 
ATOM    188 CA   ALA C  26      86.402  48.355   4.259   0.000   0.000   0.000 
ATOM    193 CA   THR C  27      88.767  45.715   5.625   3.049   0.000   0.000 
ATOM    200 CA   VAL C  28      86.811  45.515   8.869   2.787   0.000   0.000 
ATOM    207 CA   LYS C  29      83.557  45.184   6.898   4.328   0.000   0.000 
ATOM    216 CA   ALA C  30      84.997  42.387   4.777   6.413   0.000   0.000 
ATOM    221 CA   LYS C  31      86.121  40.522   7.889   6.621   0.000   0.000 
ATOM    230 CA   LYS C  32      82.787  41.128   9.641   6.766   0.000   0.000 
ATOM    239 CA   GLU C  33      81.053  39.523   6.681   1.910   0.000   0.000 
ATOM    248 CA   ARG C  34      83.361  36.501   6.568   0.000   0.000   0.000 
ATOM    259 CA   ASP C  35      82.888  36.000  10.284   0.000   0.000   0.000 
ATOM    267 CA   ALA C  36      79.146  36.539  10.043   0.000   0.000   0.000 



