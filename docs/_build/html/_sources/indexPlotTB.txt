:tocdepth: 2

PlotTB
========================

.. toctree::
    :hidden:

    codePlotTB


**plotTB** plots the **TB** outputs.

**plotTB** can:

    * plot lattice with bonds. 

    * plot spectrum with/without a selected sublattice polarization. 

    * plot states. 


Examples
--------------------------------------

Dimer chain:

.. code::

    from latticeTB import *
    from eigTB import *
    from plotTB import *
    from propTB import *

    nx, ny = 9, 1
    ri = [[0, 0], [1, 0]]
    tags = [b'a', b'b']
    lat_chain = latticeTB(tags=tags, ri=ri, nor=2, ang=0)
    lat_chain.get_lattice(nx=nx, ny=ny)
    lat_chain.remove_sites([17])
    eig_chain = eigTB(lat_chain)
    t1 = 1.
    eig_chain.set_onsite([0, -.5j])
    eig_chain.set_hop([t1])
    eig_chain.get_ham()
    eig_chain.get_eig(eigenvec=True)
    zero_mode = eig_chain.get_state_pola(pola_tag=b'a')

    plt_chain = plotTB(eig_chain)
    fig_lat = plt_chain.plt_lattice(ms=12, figsize=(8, 1))
    fig_spec = plt_chain.plt_spec(pola_tag=b'a')
    fig_zero_mode = plt_chain.plt_intensity1d(zero_mode)

    prop = propTB(lat=lat_chain, steps=800, dz=0.05)
    psi_init = np.ones(eig_chain.sites, 'c16') / np.sqrt(eig_chain.sites)
    prop.get_prop(ham=eig_chain.ham, psi_init=psi_init, norm=True)
    fig_prop = prop.plt_prop1d()
    plt.show()

    save_chain = saveFigTB(sys=eig_chain, dir_name='chain', params={'t1': t1})
    save_chain.save_fig_lat(fig_lat, 'lattice')
    save_chain.save_fig(fig_spec, 'spec')
    save_chain.save_fig(fig_zero_mode, 'zero_mode')
    save_chain.save_fig(fig_prop, 'prop')


.. image:: ../TBfig/chain_n17/lattice.png
    :height: 100px
    :width:  95%
    :align: center

.. image:: ../TBfig/chain_n17/spec_ea0j_eb-0,5j.png
    :height: 100px
    :width:  55%
    :align: center

.. image:: ../TBfig/chain_n17/zero_mode_ea0j_eb-0,5j.png
    :height: 100px
    :width:  55%
    :align: center
.. image:: ../TBfig/chain_n17/prop_ea0j_eb-0,5j.png
    :height: 100px
    :width:  65%
    :align: center

Line centrered hexagonal lattice:

.. code::

    from latticeTB import *
    from eigTB import *
    from plotTB import *
    from math import pi
    from collections import OrderedDict

    nx, ny = 5, 5
    ri = [[0, 0], [0.5*sqrt(3), 0.5], [0.25*sqrt(3), 0.25], [0.75*sqrt(3), 0.25], [0.5*sqrt(3), 1.]]
    tags = [b'a', b'b', b'c', b'd', b'e']
    hlc = latticeTB(tags=tags, ri=ri, nor=sqrt(3), ang=pi/3)
    hlc.get_lattice(nx=nx, ny=ny)
    hlc.remove_dangling(nor_bond=0.5)
    eig_hlc = eigTB(hlc)
    t1, t2 = 1., .5
    eig_hlc.set_hop([t1, t2])
    eig_hlc.get_ham()
    eig_hlc.get_eig(eigenvec=True)
    zero_mode = eig_hlc.get_state_pola(pola_tag=b'a')
    branch_neg = eig_hlc.get_states_en(e_min=-6, e_max=-2.1)
    flat_band = eig_hlc.get_states_en(e_min=-2.01, e_max=-1.99)
    branch_pos = eig_hlc.get_states_en(e_min=.1, e_max=6)

    plt_hlc = plotTB(eig_hlc)
    fig_lat_hlc = plt_hlc.plt_lattice(ms=12)
    fig_spec_hlc = plt_hlc.plt_spec(pola_tag=b'a')
    fig_zero_mode_hlc = plt_hlc.plt_intensity(zero_mode, title='')
    fig_banch_neg_hlc = plt_hlc.plt_intensity(branch_neg, title='')
    fig_flat_band_hlc = plt_hlc.plt_intensity(flat_band, title='')
    fig_banch_pos_hlc = plt_hlc.plt_intensity(branch_pos, title='')
    fig_branch_neg_disk_hlc = plt_hlc.plt_intensity_disk(branch_neg, s=2000, title='')
    plt.show()

    save_hlc = saveFigTB(sys=hlc, dir_name='hexa_lc', params=OrderedDict([('t1', t1), ('t2', t2)]))
    save_hlc.save_fig_lat(fig_lat_hlc, 'lattice')
    save_hlc.save_fig(fig_spec_hlc, 'spec')
    save_hlc.save_fig(fig_zero_mode_hlc, 'zero_mode')
    save_hlc.save_fig(fig_banch_neg_hlc, 'banch_neg')
    save_hlc.save_fig(fig_flat_band_hlc, 'flat_band')
    save_hlc.save_fig(fig_banch_pos_hlc, 'banch_pos')
    save_hlc.save_fig(fig_branch_neg_disk_hlc, 'banch_neg_disk')


.. image:: ../TBfig/hexa_lc_n111/lattice.png
    :height: 100px
    :width:  65%
    :align: center

.. image:: ../TBfig/hexa_lc_n111/spec_t1(1+0j)_t2(0,5+0j).png
    :height: 100px
    :width:  55%
    :align: center

.. image:: ../TBfig/hexa_lc_n111/zero_mode_t1(1+0j)_t2(0,5+0j).png
    :height: 100px
    :width:  55%
    :align: center

.. image:: ../TBfig/hexa_lc_n111/banch_neg_t1(1+0j)_t2(0,5+0j).png
    :height: 100px
    :width:  55%
    :align: center

.. image:: ../TBfig/hexa_lc_n111/banch_neg_disk_t1(1+0j)_t2(0,5+0j).png
    :height: 100px
    :width:  55%
    :align: center

.. image:: ../TBfig/hexa_lc_n111/flat_band_t1(1+0j)_t2(0,5+0j).png
    :height: 100px
    :width:  55%
    :align: center

.. image:: ../TBfig/hexa_lc_n111/banch_pos_t1(1+0j)_t2(0,5+0j).png
    :height: 100px
    :width:  55%
    :align: center


Feedback
^^^^^^^^^^^^^^

Please send comments or suggestions for improvement to cpoli83 at hotmail dot fr