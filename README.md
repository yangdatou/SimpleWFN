# SimpleWFN
This project closely follows [ProgrammingProjects](https://github.com/CrawfordGroup/ProgrammingProjects/tree/master/Project%2301) from Crawford Group. 
The purpose of this project is to practice fundamental C++ programming techniques in the context of electronic structure problems and to provide a deeper
understanding of Hartree-Fock theory and post-HF algorithms.

## Hartree-Fock Theory
The theoretical background can be found in Ch. 3 of the text by Szabo and Ostlund or in the 
[nice set of on-line notes](http://vergil.chemistry.gatech.edu/notes/hf-intro/hf-intro.html) written by David Sherrill.

- Step 1: Nuclear Repulsion Energy 

- Step 2: One-Electron Integrals

- Step 3: Two-Electron Integrals

- Step 4: Build the Orthogonalization Matrix

- Step 5: Build the Initial Guess Density

- Step 6: Compute the Inital SCF Energy

- Step 7: Compute the New Fock Matrix 

- Step 8: Build the New Density Matrix 
<!--
## Step 5: Build the Initial Guess Density

Form an initial (guess) Fock matrix in the orthonormal AO basis using the core Hamiltonian as a guess:



Diagonalize the Fock matrix:



Note that the &epsilon;<sub>0</sub> matrix contains the initial orbital energies.

Transform the eigenvectors into the original (non-orthogonal) AO basis:



Build the density matrix using the occupied MOs:



where *m* indexes the columns of the coefficient matrices, and the summation includes only the occupied spatial MOs.

  * [Hint 1](./hints/hint5-1.md): Transformed Fock matrix
  * [Hint 2](./hints/hint5-2.md): Initial MO Coefficients
  * [Hint 3](./hints/hint5-3.md): Initial Density Matrix

## Step 6: Compute the Inital SCF Energy

The SCF electronic energy may be computed using the density matrix as:



The total energy is the sum of the electronic energy and the nuclear repulsion energy:



where *0* denotes the initial SCF energy.

 * [Hint 1](./hints/hint6-1.md): Initial Electronic Energy

## Step #7: Compute the New Fock Matrix 

Start the SCF iterative procedure by building a new Fock matrix using the previous iteration's density as:



where the double-summation runs over all the AOs and *i-1* denotes the density for the last iteration.

  * [Hint 1](./hints/hint7-1.md): New Fock Matrix
  * [Hint 2](./hints/hint7-2.md): Fock-Build Code

## Step #8: Build the New Density Matrix 

Form the new density matrix following the same procedure as in Step #5 above:

Orthogonalize:



Diagonalize:



Back-transform:



Compute the density:



where *i* denotes the current iteration density.

## Step #9: Compute the New SCF Energy 

Compute the new SCF energy as before:



where *i* denotes the SCF energy for the *i*th iteration.

## Step #10: Test for Convergence 
Test both the energy and the density for convergence:



If the difference in consecutive SCF energy and the root-mean-squared difference in consecutive densities do not fall below the prescribed thresholds, return to Step #7 and continue from there.

  * [Hint 1](./hints/hint10-1.md): Energies for Each Iteration

## Additional Concepts
###  The MO-Basis Fock Matrix
At convergence, the canonical Hartree-Fock MOs are, by definition, eigenfunctions of the Fock operator, viz.



If we multiply on the left by an arbitrary MO and integrate, we obtain:



In other words, the Fock matrix should be diagonal in the MO basis, with the orbital energies as its diagonal elements.  We can demonstrate this explicitly using the AO-basis Fock matrix by first re-writing the above expression using the LCAO-MO coefficients:



Use the above equation to transform the Fock matrix from the AO basis to the MO basis and demonstrate that it is indeed diagonal (to within the convergence limits of the SCF iterative procedure).

### One-Electron Properties 
As discussed in detail in Ch. 3 of the text by Szabo and Ostlund, the calculation of one-electron properties requires density matrix and the relevant property integrals.  The electronic contribution to the electric-dipole moment may be computed using,



where the vector notation implies three sets of dipole-moment integrals -- one for each Cartesian component of the dipole operator.

Two points to note:
  - In order to compute the total dipole moment, you must include the nuclear contribution, which requires the atomic numbers and Cartesian coordinates of the nuclei, in addition to the above.
  - The factor 2 appearing above arises because the definition of the density used in this project differs from that used in Szabo & Ostlund.

The test cases provided below include the structural information dipole integrals needed to compute the dipole moment.

### Population Analysis/Atomic Charges
A Mulliken population analysis (also described in Szabo & Ostlund, Ch. 3) requires the overlap integrals and the electron density, in addition to information about the number of basis functions centered on each atom.  The charge on atom *A* may be computed as:



where the summation is limited to only those basis functions centered on atom *A*.
-->


## Run Code
To run the code, use the following commands:
```
git clone git@github.com:yangdatou/SimpleWFN.git
mkdir bin
make test
```

## References
Szabo and Ostlund. *Modern Quantum Chemistry: Introduction to Advanced Electronic Structure Theory*, Dover Publications (1996)