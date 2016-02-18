Python Tight-Binding module
========================


![Alt text](https://github.com/cpoli/tbee/blob/master/logoTBee.png)



**Python** module to build up and solve **Tight-Binding** models. 

**tbee** is written in fully vectorized **Numpy**.

**tbee** is composed of the following classes:

    * lattice
    * system
    * plot
    * propagation
    * save


**tbee** main features:

    * Complex lattice structures.
    * Complex-valued onsite energies and hoppings.
    * Hermitian and non-Hermitian Tight-Binding Hamiltonians.
    * Sublattices.
    * Hoppings defined by their type, tags, and angles.
    * Any type of hoppings:

        * Neighbors hoppings,
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

**tbee** is available at https://github.com/cpoli/tbee


To use **tbee**:

  * Install Python3.4 or Python3.5 and three additional packages:

      * numpy 1.10
      * scipy 0.16
      * matplotlib 1.5

  * See https://cpoli.github.io/python-doc.html for Python installation details
      and to install a github repository (for mac).

Examples are available at https://github.com/cpoli/tbee/tree/master/examples