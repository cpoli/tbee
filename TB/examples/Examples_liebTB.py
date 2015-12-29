from latticeTB import *
from plotTB import *
from eigTB import *
from liebTB import *


n1, n2 = 5, 5
s = 200
alpha = 0.3
lat = liebTB()
lat.get_lattice(n1=n1, n2=n2)
lat.remove_dangling()
sys = liebEig(lat)
plot = plotTB(sys)

t_ab, t_ba, t_ac, t_ca = np.round(4/3, 2), np.round(2/3, 2), np.round(4/3, 2), np.round(2/3, 2)
#c = 0.5
#t_pp = c * np.sqrt(t_ba + t_ac)
#t_pm = c * np.sqrt(t_ba + t_ca)
#t_mp = c * np.sqrt(t_ab + t_ac)
#t_mm = c * np.sqrt(t_ab + t_ca)
t_mp, t_pp, t_mm, t_pm = 0.4, 0.3, 0.3, 0.1

sys.set_nearest_neighbor_hop(t_ab=t_ab, t_ba=t_ba, t_ac=t_ac, t_ca=t_ca)
#sys.set_ons({b'a': 0j, b'b': -0.2j, b'c': -0.2j})
sys.set_next_nearest_neighbor_hop(t_pp=t_pp, t_pm=t_pm, t_mp=t_mp ,t_mm=t_mm)
sys.set_defect_dimer_x(4)
sys.set_defect_dimer_y(4)
#sys.set_disorder_placket(alpha)
#sys.set_disorder_generic(alpha)
#sys.set_disorder_pair(alpha)
sys.get_ham(complex_transpose=True)
sys.get_eig(eigenvec=True)

sys.get_coor_hop_dis()


zero_mode = sys.get_intensity_pola(tag_pola=b'a')

fig_spectrum = plot.spectrum(tag_pola=b'a', fs=25, ms=10)
fig_zero_mode = plot.intensity_disk(zero_mode, s=s, fs=25)
fir_lat_hop = plot.lattice_hop(ms=15)


'''
from collections import OrderedDict
save = saveTB(sys=sys, params=OrderedDict([('ta', t_ab), ('tb', t_ba), ('tc', t_ca), ('td', t_ca)]), 
          dir_name='lieb_supp', ext='pdf')
save.fig(fig_zero_mode, 'zero_mode_uni')
save.fig(fig_spectrum, 'spectrum_uni')
#save.fig(fig_coor_real, 'lattice_realXY')


#sys.hop['t'] = 2*np.exp(-sys.hop['t'])
if alpha == 0:
    sys.get_coor_hop()
else:
    sys.get_coor_hop_dis()
fig_coor_real = plot.lattice_hop()

from collections import OrderedDict
save = saveTB(sys=sys, params=OrderedDict([('ta', t_ab), ('tb', t_ba), ('tc', t_ca), ('td', t_ca)]), 
                                     dir_name='lieb_supp', ext='pdf')
save.fig(fig_coor_real, 'lattice_real_gen2')


sys.rename_hop_tag([{'n':1, 'ang':90, 'tag_new': b'u'}, 
                                    {'n':1,'ang':30, 'tag_new': b'v'}, 
                                    {'n':1,'ang':150, 'tag_new': b'w'}])
'''
plt.show() 