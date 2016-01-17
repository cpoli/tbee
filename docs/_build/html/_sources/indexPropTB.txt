:tocdepth: 2

Prop TB
========================

.. toctree::
    :hidden:

    codePropTB


**propTB** gets the field Time Evolution.

**propTB** can:

    * consider linear time dependent Schrödinger equation. 

    * consider non-linear time dependent Schrödinger equation. 


Examples
--------------------------------------

Lieb Lattice - Passive amplication of a zero-mode:

.. code::

    from latticeTB import *
    from eigTB import *
    from plotTB import *
    from propTB import *
    from math import pi
    from collections import OrderedDict
    
    nx, ny = 9, 9
    ri = [[0, 0], [1, 0], [0, 1]]
    tags = [b'a', b'b', b'c']
    lat_lieb = latticeTB(tags=tags, ri=ri, nor=2, ang=pi/2)
    lat_lieb.get_lattice(nx=nx, ny=ny)
    lat_lieb.remove_dangling(nor_bond=1.)

    eig_lieb = eigTB(lat_lieb)
    t1, t2 = 1., 2.
    eig_lieb.set_hop([t1, t2])
    eig_lieb.set_onsite([0, -.2j, -.2j])
    eig_lieb.get_ham()
    eig_lieb.get_eig(eigenvec=True)
    zero_mode = eig_lieb.get_state_pola(pola_tag=b'a')

    plt_lieb = plotTB(eig_lieb)
    fig_lat = plt_lieb.plt_lattice(ms=15)
    fig_spec = plt_lieb.plt_spec(pola_tag=b'a')
    fig_zero_mode = plt_lieb.plt_intensity(zero_mode)

    prop = propTB(lat=lat_lieb, steps=150, dz=0.05)
    psi_init = np.ones(eig_lieb.sites, 'c16') / np.sqrt(eig_lieb.sites)
    prop.get_prop(ham=eig_lieb.ham, psi_init=psi_init, norm=True)
    ani = prop.get_ani(s=200)
    plt.show()

    save_lieb = saveFigTB(sys=eig_lieb, dir_name='lieb', 
                                          params=OrderedDict([('t1', t1), ('t2', t2)]))
    save_lieb.save_fig_lat(fig_lat, 'lat')
    save_lieb.save_fig_lat(fig_spec, 'spec')
    save_lieb.save_fig(fig_zero_mode, 'zero_mode')
    save_lieb.save_ani(ani, 'ani')


.. image:: ../TBfig/lieb_n225/lat.png
    :height: 100px
    :width:  45%

.. image:: ../TBfig/lieb_n225/spec_ea0j_eb-2j_ec-2j_t1(1+0j)_t2(0,2+0j).png
    :height: 100px
    :width:  45%

.. image:: ../TBfig/lieb_n225/zero_mode_ea0j_eb-2j_ec-2j_t1(1+0j)_t2(0,2+0j).png
    :height: 100px
    :width:  55%
    :align: center


Feedback
^^^^^^^^^^^^^^

Please send comments or suggestions for improvement to cpoli83 at hotmail dot fr