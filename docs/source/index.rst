.. tbee documentation master file, created by
   sphinx-quickstart on Sat Feb 13 09:54:24 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

tbee
====

.. toctree::
    :hidden:

    tutorial
    history
    api/gallery/index
    tbee

**tbee** is a Python package to build and solve **Tight-Binding** models,
written in fully vectorized **NumPy**. It's built for physics students
working through a textbook problem, curious learners exploring a topic
on their own, and educators building a demonstration.

Start with the :doc:`tutorial`, browse the :doc:`example gallery
</api/gallery/index>`, or jump straight to the :doc:`tbee` API reference.

**tbee** is composed of the following classes and modules:

    * Lattice
    * System
    * KSpace
    * Plot
    * Propagation
    * Save
    * lattices
    * dos


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
       * Magnetic field (Peierls substitution).
       * Strain.
       * Hopping disorder.
       * Onsite disorder.

    * Reciprocal-space Bloch Hamiltonians and band structures.
    * Berry curvature and Chern numbers.
    * An optional spin-1/2 degree of freedom, for spin-orbit coupling and
      Zeeman terms.
    * Ribbons (edge states) cut from any periodic model.
    * Broadened density of states.
    * A small library of ready-made lattices.
    * Time propagation.

**tbee** is available at https://github.com/cpoli/tbee


To use **tbee**:

  * Install Python 3.10 or later and three additional packages:

      * numpy
      * scipy
      * matplotlib

  * ``pip install -e .`` from a clone of the repository.

Examples are available at https://github.com/cpoli/tbee/tree/master/examples,
and rendered with their output in the :doc:`example gallery
</api/gallery/index>`.