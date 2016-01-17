:tocdepth: 2

Tight-Binding
==============================

.. toctree::
    :hidden:

    indexLattice
    indexSystem
    indexPlot
    indexPropagation


.. image:: logoTB.png
    :width:  35%
    :scale: 50 %
    :align: center

**TB** is a **Python** module building and solving **Tight-Binding** models. 

Sphinx documentation of **TB**: http://cpoli.github.io/TBdoc/indexTB.html

**TB** is composed of:

    * lattice
    * system
    * plot
    * propagation
    * save

**TB** is written in fully vectorized **Numpy**.:

**TB** main features:

    * Complex lattice structures.
    * Complex-valued onsite energies and hoppings.
    * Hermitian and non-Hermitian Tight-Binding Hamiltonians.
    * Onsite energies defined by tags.
    * Any type of hoppings:
        * Neighbors hoppings
        * Next-neighbors hoppings, 
        * Next-next-neighbors hoppings,
        * etc..
    * Hoppings defined by their type, tags, and angles.
    * Implementation of onsite energies and hopping patterns:
       * Dimerization defects (change of hopping patterns).
       * Implementation of magnetic field.
       * Implementation of strain.
       * Local value changes. 
       * Hopping disorder.
       * Onsite disorder.
    * Time propagation
   
Installation
----------------

**TB** is available at https://github.com/cpoli/TB

To use **TB**,  you need to install the programming language python and three additional packages:

    * python 3.x
    * numpy
    * scipy
    * matplotlib

See https://cpoli.github.io/python-doc.html for Python installation details.

Feedback
-----------------

Please send comments or suggestions for improvement to cpoli83 at hotmail dot fr