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

## The Second-Order Moller-Plesset Perturbation (MP2)
Unlike the Hartree-Fock energy, corrlation methods like the MP2 energy are usually expressed in terms of MO-basis quantities (integrals, MO energies).
One of the expensive part of the calculation is the transformation of the two-electron integrals from the AO basis to MO basis.
The purpose of this project is to understand this transformation and fundamental techniques for its efficient implementation.

- Step 1: Transform the two-electron integrals to MO basis: (a) the noddy algorithm, and (b) the smarter algorithm.
- Step 2: Compute the MP2 amplitude and MP2 energy

## Coupled-Cluster Singles and Doubles (CCSD) 
The coupled cluster model provides a higher level of accuracy beyond the MP2 approach. The purpose of this project is to understand the fundamental aspects of the calculation of the CCSD wave function and corresponding energy. 
## Run Code
To run the code, use the following commands:
```
git clone git@github.com:yangdatou/SimpleWFN.git
mkdir bin
make test
```

## References
- Szabo and Ostlund. *Modern Quantum Chemistry: Introduction to Advanced Electronic Structure Theory*, Dover Publications (1996)
- J. Chem. Phys. 94, 4334 (1991); https://doi.org/10.1063/1.460620
- J. Chem. Phys. 120, 2581 (2004); https://doi.org/10.1063/1.1637577
- PySCF: https://github.com/pyscf/pyscf