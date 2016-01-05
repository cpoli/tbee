from lattice import *
from system import *
from plot import *

unit_cell = [{'tag': b'a', 'r0': [0, 0]}, 
             {'tag': b'b', 'r0': [1, 0]}]
prim_vec = {'norm': 2, 'angle': 0}

lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
sys = system(lat)
plt = plot(sys=sys)
n1, n2 = 10, 1
lat.get_lattice(n1=n1, n2=n2)

a = np.arange(20)
i = a % 2 == 0
b = np.delete(a, np.where(a % 2 == 0))

sys.set_hopping([{'n': 1, 't': 1.}])
sys.get_ham()
sys.get_eig(eigenvec=True)
sys.get_ipr()

'''
class A:
    def __init__(self, a):
        self.a = a
        self.v = np.array([9])

    def vv(self):
        self.v = np.array([999])
class B:
    def __init__(self, a, b):
        self.a = a
        self.b = b

aa = A(1)
bb = B(aa, 2)
print(bb.a.a)
aa.a = 3
print(bb.a.a)
print(bb.a.v)
aa.vv()
print(aa.v)
print(bb.a.v)
'''