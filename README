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
 Kumar, S. and Bansal, M. (1996). Structural and Sequence Characteristics of Long Alpha Helices in Globular Proteins. Biophysical J.,71, 1574-1586.

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

    biHelix [-i=required] [-o=helix_out.txt] [-k] [-q] [-v]  

Options and flags {default values}:
    -i, --in    The input file name.  {required}  
    -o, --out    The output file name.  {helix_out.txt}
    -k, --kink    Only calculate the kinks.
    -q, --quiet    Only output to file.
    -v, --verbose    write more to file.

  Positional arguments: None

Also, -?, -h, -H, -help, --help, and --usage are recognized.


Input format
-------------
 see ATOM section used in PDB format

Output format
--------------
 see ATOM section used in PDB format.
