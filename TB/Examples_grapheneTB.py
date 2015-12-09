from grapheneTB import *
from liebTB import *

'''
# Quick example
lat = grapheneTB()
lat.get_lattice(n1=15, n2=15)
lat.remove_dangling()

sys = grapheneEig(lat=lat)
plot = plotTB(sys=sys)
print(lat.coor)
sys.print_hop(n=2)
plot.lattice(ms=8);
sys.set_hop_uni([{'n':1, 't':1.}, {'n':2, 't': 0.001}])
sys.print_hop(n=2)

sys.get_ham()
sys.get_eig()
plot.spectrum()
'''
'''
# Strained Hexagonal Flakes
lat = grapheneTB()
lat.hexagon_zigzag(5)
sys = grapheneEig(lat)
plot = plotTB(sys)
betas = sys.get_beta_lims()

# increasing negative strain 
sys.set_hop_linear_strain(t=1, beta=0)
sys.get_coor_hop()
plot.lattice_hop(ms=5, lw=1, plt_hop=True, figsize=(3, 3));

sys.set_hop_linear_strain(t=1, beta=1/4*betas[0])
sys.get_coor_hop()
plot.lattice_hop(ms=5, lw=1, plt_hop=True, figsize=(3, 3));

sys.set_hop_linear_strain(t=1, beta=2/4*betas[0]);
sys.get_coor_hop()
plot.lattice_hop(ms=5, lw=1, plt_hop=True, figsize=(3, 3));

sys.set_hop_linear_strain(t=1, beta=3/4*betas[0])
sys.get_coor_hop()
plot.lattice_hop(ms=5, lw=1, plt_hop=True, figsize=(3, 3));

# increasing positive strain 
sys.set_hop_linear_strain(t=1, beta=0)
sys.get_coor_hop()
plot.lattice_hop(ms=5, lw=1, plt_hop=True, figsize=(3, 3));

sys.set_hop_linear_strain(t=1, beta=1/4*betas[1])
sys.get_coor_hop()
plot.lattice_hop(ms=5, lw=1, plt_hop=True, figsize=(3, 3));

sys.set_hop_linear_strain(t=1, beta=2/4*betas[1])
sys.get_coor_hop()
plot.lattice_hop(ms=5, lw=1, plt_hop=True, figsize=(3, 3));

sys.set_hop_linear_strain(t=1, beta=3/4*betas[1])
sys.get_coor_hop()
plot.lattice_hop(ms=5, lw=1, plt_hop=True, figsize=(3, 3));

# Strained Triangular Flakes
lat = grapheneTB()
lat.triangle_zigzag(10)
sys = grapheneEig(lat)
plot = plotTB(sys)
betas = sys.get_beta_lims()

# increasing negative strain 
sys.set_hop_linear_strain(t=1, beta=0)
sys.get_coor_hop()
plot.lattice_hop(ms=5, lw=1, plt_hop=True, figsize=(3, 3));

sys.set_hop_linear_strain(t=1, beta=1/4*betas[0])
sys.get_coor_hop()
plot.lattice_hop(ms=5, lw=1, plt_hop=True, figsize=(3, 3));

sys.set_hop_linear_strain(t=1, beta=2/4*betas[0])
sys.get_coor_hop()
plot.lattice_hop(ms=5, lw=1, plt_hop=True, figsize=(3, 3));

sys.set_hop_linear_strain(t=1, beta=3/4*betas[0])
sys.get_coor_hop()
plot.lattice_hop(ms=5, lw=1, plt_hop=True, figsize=(3, 3));

# increasing positive strain 
sys.set_hop_linear_strain(t=1, beta=0)
sys.get_coor_hop()
plot.lattice_hop(ms=5, lw=1, plt_hop=True, figsize=(3, 3));

sys.set_hop_linear_strain(t=1, beta=1/4*betas[1])
sys.get_coor_hop()
plot.lattice_hop(ms=5, lw=1, plt_hop=True, figsize=(3, 3));

sys.set_hop_linear_strain(t=1, beta=2/4*betas[1])
sys.get_coor_hop()
plot.lattice_hop(ms=5, lw=1, plt_hop=True, figsize=(3, 3));

sys.set_hop_linear_strain(t=1, beta=3/4*betas[1])
sys.get_coor_hop()
plot.lattice_hop(ms=5, lw=1, plt_hop=True, figsize=(3, 3));

# PSEUDO LANDAU LEVELS

# Triangular flake with zigzag terminations
lat = grapheneTB()
lat.triangle_zigzag(24)
sys = grapheneEig(lat)
plot = plotTB(sys)
betas = sys.get_beta_lims()
sys.set_hop_linear_strain(t=1, beta=betas[1])
sys.get_coor_hop()
plot.lattice(ms=5);
plot.lattice_hop(ms=5);

sys.get_ham(complex_transpose=True)
sys.get_eig(eigenvec=True)
plot.spectrum(tag_pola=b'a', en_lims=[-2.2, 2.2]);

ll0 = sys.get_intensity_en(en_lims=[-0.1, 0.1])
ll1 = sys.get_intensity_en(en_lims=[0.8, 1.2])
ll2 = sys.get_intensity_en(en_lims=[1.2, 1.4])
ll3 = sys.get_intensity_en(en_lims=[1.4, 1.6])
ll4 = sys.get_intensity_en(en_lims=[1.6, 1.8])
ll5 = sys.get_intensity_en(en_lims=[1.8, 2.])
plot.intensity_area(ll0, s=20, plt_hop=True);
plot.intensity_area(ll1, s=20, plt_hop=True);
plot.intensity_area(ll2, s=20, plt_hop=True);




lat = grapheneTB()
lat.triangle_zigzag(24)
lat.remove_sites([lat.sites//2])
sys = grapheneEig(lat)
plot = plotTB(sys)
betas = sys.get_beta_lims()
sys.set_hop_linear_strain(t=1, beta=betas[1])
plot.lattice(ms=5);
sys.get_ham(complex_transpose=True)
sys.get_eig(eigenvec=True)
plot.spectrum(ms=8, tag_pola=b'a', en_lims=[-2.2, 2.2]);
plot.show()
'''


#unit_cell = [{'tag': b'a', 'r0': [0, 0]}]
unit_cell = [{'tag': b'a', 'r0': [0, 0]}, 
                   {'tag': b'b', 'r0': [1, 0]},
                   {'tag': b'c', 'r0': [0, 1]}]
prim_vec = {'norm': 2, 'angle': 90}
lat = latticeTB(unit_cell=unit_cell, prim_vec=prim_vec)
lat.get_lattice(n1=5, n2=5)
lat.remove_dangling()

sys = eigTB(lat=lat)
plot = plotTB(sys=sys)
sys.print_hop(3)

#sys.set_hopping([{'n': 1}, {'n': 1}, {'n': 1}, {'n': 4}, {'n': 4}, {'n': 5}, {'n': 2}])
#sys.set_hopping([{'n': 1, 'ang': 0, 't':2}, 
#                            {'n': 1, 'ang': 90, 't':3}, 
#                           {'n': 2, 'ang': 45, 't': 2}])
# nearest hopping
sys.set_hopping([{'n': 1, 'tag': b'ab', 't':2}, 
                            {'n': 1, 'tag': b'ba', 't':1}, 
                            {'n': 1, 'tag': b'ac', 't':2},
                            {'n': 1, 'tag': b'ca', 't':1}])
# next nearest
sys.set_hopping([{'n': 2, 'ang': 45, 'tag': b'bc', 't': 1.}, 
                            {'n': 2, 'ang': 135, 'tag': b'bc', 't': 1.5}, 
                            {'n': 2, 'ang': 45, 'tag': b'cb', 't': 1.},
                            {'n': 2, 'ang': 135, 'tag': b'cb', 't': .5}])

#sys.set_defect_dimer_y(5)
sys.get_ham(complex_transpose=True)
sys.get_eig(eigenvec=True)
sys.get_intensity_pola(b'a')
plot.spectrum()

plot.lattice(ms=10, plt_hop=True)
plot.show()
