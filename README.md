Python Tight-Binding module
========================


![Alt text](https://github.com/cpoli/TB/blob/master/logoTB_.png)



**Python** module to build up and solve **Tight-Binding** models. 

**TB** is written in fully vectorized **Numpy**

**TB** is composed of the following classes:

    * lattice
    * system
    * plot
    * propagation
    * save


**TB** main features:

    * Complex lattice structures.
    * Complex-valued onsite energies and hoppings.
    * Hermitian and non-Hermitian Tight-Binding Hamiltonians.
    * Onsite energies defined by tags.
    * Hoppings defined by their type, tags, and angles.
    * Any type of hoppings:
        * Neighbors hoppings
        * Next-neighbors hoppings, 
        * Next-next-neighbors hoppings,
        * etc..
    * Implementation of onsite energies and hopping patterns:
       * Dimerization defects.
       * Magnetic field.
       * Strain.
       * Hopping disorder.
       * Onsite disorder.
    * Time propagation

**TB** is available at https://github.com/cpoli/TB


To use TB:

  * Install Python3.5 and three additional packages:

      * numpy 1.10
      * scipy 0.16
      * matplotlib 1.5

    * See https://cpoli.github.io/python-doc.html for Python installation details
      and to install a github repository (for mac).

Examples are available at https://github.com/cpoli/TB/tree/master/examples