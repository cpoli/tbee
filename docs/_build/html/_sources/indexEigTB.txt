:tocdepth: 2

Eig TB
========================

.. toctree::
    :hidden:

    codeEigTB


**eigTB** solves the Tight-Binding Hamiltonian.

**eigTB** can:

    * set the lattice hoppings. 

    * set the lattice onsite energies. 

    * set hopping and onsite disorder. 

    * solve the Tight-Binding Hamiltonian:

        * get eigenenergies and eigenvectors.

        * get sublattice polarizations.

    * select states using energy range or polarization condition.


Examples
--------------------------------------

Square lattice

Lattice creation:

.. code::

    # imports
    from latticeTB import *
    from eigTB import *
    from plotTB import *
    from math import pi
    from collections import OrderedDict

    # lattice
    nx, ny = 5, 5
    ri = [[0, 0]]
    tags = [b'a']
    lat = latticeTB(tags=tags, ri=ri, nor=1, ang=pi/2)
    lat.get_lattice(nx=nx, ny=ny)
    fig_lat = lat.plt_lattice(ms=18, plt_label=True)


.. image:: ../TBfig/square_n25/lat.png
    :height: 100px
    :width:  45%
    :align: center


Tight-Binding eigenvalue problem:

    *  with nearest-neighbor hoppings

        .. code::

            t1 = 1.
            eig = eigTB(lat)
            eig.set_hop(ho=[t1])
            eig.get_ham()
            eig.get_eig()
            plot = plotTB(eig) 
            fig_lat_bonds = plot.plt_lattice(ms=18)
            fig_spec = plot.plt_spec()
            plt.show()
            # save 
            save = saveFigTB(sys=eig, dir_name='square', params={'t1': t1})
            save.save_fig_lat(fig_lat, 'lat')
            save.save_fig(fig_lat_bonds, 'lattice')
            save.save_fig(fig_spec, 'spec')


        .. image:: ../TBfig/square_n25/lattice_t1(1+0j).png
            :height: 100px
            :width:  45%

        .. image:: ../TBfig/square_n25/spec_t1(1+0j).png
            :height: 100px
            :width:  45%

    * adding next nearest-neighbor hoppings

        .. code::

            t2 = .8
            eig.set_hop(ho=[t1, t2])
            eig.get_ham()
            eig.get_eig(eigenvec=True)
            plot = plotTB(eig)
            zero = eig.get_states_en(-.001, .001)
            fig_lat_bonds = plot.plt_lattice(ms=18)
            fig_spec = plot.plt_spec()
            fig_zero = plot.plt_intensity(zero)
            plt.show()
            # save 
            save = saveFigTB(sys=eig, dir_name='square', params=OrderedDict([('t1', t1), ('t2', t2)]))
            save.save_fig(fig_lat_bonds, 'lattice')
            save.save_fig(fig_spec, 'spec')
            save.save_fig(fig_zero, 'zero_mode')


        .. image:: ../TBfig/square_n25/lattice_t1(1+0j)_t2(0,8+0j).png
            :height: 100px
            :width:  45%

        .. image:: ../TBfig/square_n25/spec_t1(1+0j)_t2(0,8+0j).png
            :height: 100px
            :width:  45%

        The next-neighbors hoppings lead to a mode of zero energy. 

        .. image:: ../TBfig/square_n25/zero_mode_t1(1+0j)_t2(0,8+0j).png
            :height: 100px
            :width:  55%
            :align: center

    * adding next next nearest-neighbor hoppings

        .. code::

            t3 = .6
            eig.set_hop(ho=[t1, t2, t3])
            eig.get_ham()
            eig.get_eig()
            plot = plotTB(eig) 
            fig_lat_bonds = plot.plt_lattice(ms=18)
            fig_spec = plot.plt_spec()
            plt.show()
            # save 
            save = saveFigTB(sys=eig, dir_name='square', params=OrderedDict([('t1', t1), ('t2', t2), ('t3', t3)]))
            save.save_fig(fig_lat_bonds, 'lattice')
            save.save_fig(fig_spec, 'spec')

        .. image:: ../TBfig/square_n25/lattice_t1(1+0j)_t2(0,8+0j)_t3(0,6+0j).png
            :height: 100px
            :width:  45%

        .. image:: ../TBfig/square_n25/spec_t1(1+0j)_t2(0,8+0j)_t3(0,6+0j).png
            :height: 100px
            :width:  45%

    * adding next next next nearest-neighbor hoppings

        .. code::

            t4 = .4
            eig.set_hop(ho=[t1, t2, t3, t4])
            eig.get_ham()
            eig.get_eig()
            plot = plotTB(eig) 
            fig_lat_bonds = plot.plt_lattice(ms=18)
            fig_spec = plot.plt_spec()
            plt.show()
            # save 
            save = saveFigTB(sys=eig, dir_name='square', params=OrderedDict([('t1', t1), ('t2', t2), ('t3', t3), ('t4', t4)]))
            save.save_fig(fig_lat_bonds, 'lattice')
            save.save_fig(fig_spec, 'spec')


        .. image:: ../TBfig/square_n25/lattice_t1(1+0j)_t2(0,8+0j)_t3(0,6+0j)_t4(0,4+0j).png
            :height: 100px
            :width:  45%

        .. image:: ../TBfig/square_n25/spec_t1(1+0j)_t2(0,8+0j)_t3(0,6+0j)_t4(0,4+0j).png
            :height: 100px
            :width:  45%

Feedback
^^^^^^^^^^^^^^

Please send comments or suggestions for improvement to cpoli83 at hotmail dot fr