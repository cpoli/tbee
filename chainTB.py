from latticeTB import latticeTB
from eigTB import eigTB
from plotTB import plotTB
import matplotlib.pyplot as plt

'''
chain = latticeTB(ri=[[0., 0.], [1., 0.]], tags=[b'a', b'b'], nor=2., ang=0.)
chain.get_lattice(nx=10, ny=1)
chain.remove_sites([chain.sites-1])
chain.plt_lattice(ms=10, figsize=(10, 2))
sys = eigTB(chain)
sys.set_hop_uni({1: 1})
sys.get_ham()
sys.get_eig(eigenvec=True)
zero_mode = sys.get_state_pola(tag_pola=b'a')
plot = plotTB(sys)
plot.plt_lattice(ms=10, figsize=(10, 2), plt_hop=True)
plot.plt_spec(tag_pola=b'a')
plot.plt_intensity1d(zero_mode)
plot = plotTB(sys)
plot.plt_lattice(ms=10, figsize=(10, 2), plt_hop=True)
plot.plt_spec(tag_pola=b'a')
plot.plt_intensity1d(zero_mode)
'''

chain = latticeTB(ri=[[0., 0.], [1., 0.]], tags=[b'a', b'b'], nor=2., ang=0.)
chain.get_lattice(nx=10, ny=1)
chain.remove_sites([chain.sites-1])
chain.plt_lattice(ms=10, figsize=(10, 2))
sys = eigTB(chain)
sys.set_hop_nearest(dict_hop={b'ab': 2, b'ba': 1})
sys.get_ham(complex_transpose=True)
sys.get_eig(eigenvec=True)
zero_mode = sys.get_state_pola(tag_pola=b'a')
plot = plotTB(sys)
plot.plt_lattice(ms=10, figsize=(10, 2), plt_hop=True)
plot.plt_spec(tag_pola=b'a')
plot.plt_intensity1d(zero_mode)
plt.show()


