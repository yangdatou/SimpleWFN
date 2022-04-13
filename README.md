# Molecular-Geometry-Analysis
This project closely follows [ProgrammingProjects](https://github.com/CrawfordGroup/ProgrammingProjects/tree/master/Project%2301) from Crawford Group. The purpose of this project is to practice fundamental C++ programming techniques in the context of a scientific problem, viz. the calculation of the internal coordinates (bond lengths, bond angles, dihedral angles), moments of inertia, and rotational constants of a polyatomic molecule.

### Task 1: Read the Coordinate Data from Input
The input to the program is the set of Cartesian coordinates of the atoms (in bohr) and their associated atomic numbers:
- The first line above is the number of atoms (an integer), while the remaining lines contain the z-values and x-, y-, and z-coordinates of each atom (one integer followed by three double-precision floating-point numbers). 
- This [input file](./input/acetaldehyde.dat) (`acetaldehyde.dat`) along with a few other test cases can be found in this repository in the [input directory](./input).

### Task 2: Bond Lengths
Calculate the interatomic distances.

### Task 3: Bond Angles
Calculate all possible bond angles. 

### Task 4: Out-of-Plane Angles
Calculate all possible out-of-plane angles. 

### Task 5: Torsion/Dihedral Angles
Calculate all possible torsional angles. 

### Task 6: Center-of-Mass Translation
Find the center of mass of the molecule.

### Task 7: Principal Moments of Inertia
Calculate elements of the [moment of inertia tensor](http://en.wikipedia.org/wiki/Moment_of_inertia_tensor).

### Task 8: Rotational Constants
Compute the rotational constants in cm<sup>-1</sup> and MHz.

## Run Code
To run the code, use the following command:
```
make test
```

## References
E.B. Wilson, J.C. Decius, and P.C. Cross, __Molecular Vibrations__, McGraw-Hill, 1955.# molecular-geometry-analysis
