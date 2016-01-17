:tocdepth: 2

Lattice TB
========================

.. toctree::
    :hidden:

    codeLatticeTB


**latticeTB** creates finite two dimensional lattices.

**latticeTB** can:

    * create a generic lattice defined by:

        * position of the sites within the unit cell.

        * number of cells along :math:`x` and :math:`y`.

        * norm of the two primitive vectors (chosen to be egal).

        * angle between the two primitive vectors.

    * eliminate dangling sites (sites conneted to just another site).

    * eliminate selected sites.


Examples
--------------------------------------

Dimer chain

.. code::

    from latticeTB import *
    from plotTB import *

    nx, ny = 9, 1
    ri = [[0, 0], [1, 0]]
    tags = [b'a', b'b']
    lat_chain = latticeTB(tags=tags, ri=ri, nor=2, ang=0)
    lat_chain.get_lattice(nx=nx, ny=ny)
    fig_chain = lat_chain.plt_lattice(ms=10, figsize=(8, 1))
    plt.show()
    save_chain = saveFigTB(sys=lat_chain, dir_name='chain')
    save_chain.save_fig_lat(fig_chain, 'lattice')


.. image:: ../TBfig/chain_n18/lattice.png
    :height: 100px
    :width:  95%
    :align: center


Face centered square

.. code::

    from latticeTB import *
    from plotTB import *
    from math import pi

    nx, ny = 10, 10
    ri = [[0, 0], [1, 1]]
    tags = [b'a', b'b']
    lat_fcs = latticeTB(tags=tags, ri=ri, nor=2, ang=pi/2)
    lat_fcs.get_lattice(nx=nx, ny=ny)
    fig_fcs = lat_fcs.plt_lattice(ms=10)
    plt.show()
    save_fcs = saveFigTB(sys=lat_fcs, dir_name='fc_square')
    save_fcs.save_fig_lat(fig_fcs, 'lattice')

.. image:: ../TBfig/square_fc_n200/lattice.png
    :height: 100px
    :width:  55%
    :align: center

Lieb lattice

.. code::

    from latticeTB import *
    from plotTB import *
    from math import pi

    nx, ny = 6, 6
    ri = [[0, 0], [1, 0], [0, 1]]
    tags = [b'a', b'b', b'c']
    lat_lieb = latticeTB(tags=tags, ri=ri, nor=2, ang=pi/2)
    lat_lieb.get_lattice(nx=nx, ny=ny)
    fig_lieb = lat_lieb.plt_lattice(ms=10)
    lat_lieb.remove_dangling(nor_bond=1.)
    fig_lieb_dang = lat_lieb.plt_lattice(ms=10)
    plt.show()
    save_lieb = saveFigTB(sys=fig_lieb_dang, dir_name='lieb')
    save_lieb.save_fig_lat(fig_lieb, 'lattice')
    save_lieb.save_fig_lat(fig_lieb_dang, 'lattice_dang')

.. image:: ../TBfig/lieb_n96/lattice.png
    :height: 100px
    :width:  45%

.. image:: ../TBfig/lieb_n96/lattice_dang.png
    :height: 100px
    :width:  45%

Graphene

.. code::

    from latticeTB import *
    from plotTB import *
    from math import pi

    nx, ny = 7, 7
    ri = [[0, 0], [0.5*sqrt(3), 0.5]]
    tags = [b'a', b'b']
    lat_graphene = latticeTB(tags=tags, ri=ri, nor=sqrt(3), ang=pi/3)
    lat_graphene.get_lattice(nx=nx, ny=ny)
    lat_graphene.remove_dangling(nor_bond=1.)
    fig_graphene = lat_graphene.plt_lattice(ms=10, plt_label=True)
    lat_graphene.remove_sites(44)
    fig_graphene_rem = lat_graphene.plt_lattice(ms=10)
    plt.show()
    save_graphene = saveFigTB(sys=lat_graphene, dir_name='graphene')
    save_graphene.save_fig_lat(fig_graphene, 'lattice')
    save_graphene.save_fig_lat(fig_graphene_rem, 'lattice_rem')


.. image:: ../TBfig/graphene_n95/lattice.png
    :height: 100px
    :width:  45%

.. image:: ../TBfig/graphene_n95/lattice_rem.png
    :height: 100px
    :width:  45%

Kagome

.. code::

    from latticeTB import *
    from plotTB import *
    from math import pi

    nx, ny = 7, 7
    ri = [[0, 0], [1, 0], [1/2., sqrt(3)/2]]
    tags = [b'a', b'b', b'c']
    lat_kagome = latticeTB(tags=tags, ri=ri, nor=2, ang=pi/3)
    lat_kagome.get_lattice(nx=nx, ny=ny)
    fig_kagome = lat_kagome.plt_lattice(ms=10)
    plt.show()
    save_kagome = saveFigTB(sys=lat_kagome, dir_name='kagome')
    save_kagome.save_fig_lat(fig_kagome, 'lattice')

.. image:: ../TBfig/kagome_n147/lattice.png
    :height: 100px
    :width:  65%
    :align: center


Hexagon line centered

.. code::

    from latticeTB import *
    from plotTB import *
    from math import pi

    nx, ny = 5, 5
    ri = [[0, 0], [0.5*sqrt(3), 0.5], [0.25*sqrt(3), 0.25], [0.75*sqrt(3), 0.25], [0.5*sqrt(3), 1.]]
    tags = [b'a', b'b', b'c', b'd', b'e']
    hexa_lc = latticeTB(tags=tags, ri=ri, nor=sqrt(3), ang=pi/3)
    hexa_lc.get_lattice(nx=nx, ny=ny)
    fig_hexa_lc = hexa_lc.plt_lattice(colors=['b', 'b', 'r', 'r', 'r'], ms=10)
    hexa_lc.remove_dangling(nor_bond=0.5)
    fig_hexa_lc_dang = hexa_lc.plt_lattice(colors=['b', 'b', 'r', 'r', 'r'], ms=10)
    plt.show()
    save_hlc = saveFigTB(sys=hexa_lc, dir_name='hexa_lc')
    save_hlc.save_fig_lat(fig_hexa_lc, 'lat')
    save_hlc.save_fig_lat(fig_hexa_lc_dang, 'lat_dang')


.. image:: ../TBfig/hexa_lc_n111/lat.png
    :height: 100px
    :width:  45%

.. image:: ../TBfig/hexa_lc_n111/lat_dang.png
    :height: 100px
    :width:  45%


Feedback
^^^^^^^^^^^^^^

Please send comments or suggestions for improvement to cpoli83 at hotmail dot fr