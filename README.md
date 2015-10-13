# Python Tight-Binding module

**TB** is a **Python** module of **Tight-Binding** models. **TB**  builds up and solves finite tight-binding models with complex-valued onsite energies and complex-valued hoppings. 

**TB** is composed of:

* latticeTB
* eigTB
* plotTB
* propTB

Main features of **TB**:

* Written in fully vectorized **Numpy**.
* Easy lattice shape modification.
* Easy implementation of next neighbors hoppings, next next neighbors hoppings, etc..
* Time evolution (linear and non-linear)
* Easy implementation of magnetic field.
* Easy implementation of strain.
* Easy implementation of defects:
   * hopping disorder.
   * onsite disorder.
   * vacancy defects.
   * impurity defects (local change of hopping values).   
   * dimerization defects (change of hopping patterns).
