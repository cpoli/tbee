import numpy as np
import inspect


###############################
# GENERIC EXCEPTION HANDLING
###############################


def boolean(var, var_name):
    '''
    Check if *var* is a boolean.

    :raises TypeError: Parameter *var* must be a bool.
    '''
    if not isinstance(var, bool):
        raise TypeError('\n\nParameter {} must be a bool.\n'.format(var_name))


def positive_int(var, var_name):
    '''
    Check if *var* is a positive integer.

    :raises TypeError: Parameter *var* must be an integer.
    :raises ValueError: Parameter *var* must be a positive integer.
    '''
    if not isinstance(var, int):
        raise TypeError('\n\nParameter {} must be an integer.\n'.format(var_name))
    if var < 1:
        raise ValueError('\n\nParameter {} must be a positive integer.\n'.format(var_name))


def positive_int_lim(var, var_name, nmax):
    '''
    Check if *var* is a positive integer smaller than nmax.

    :raises TypeError: Parameter *var* must be an integer.
    :raises ValueError: Parameter *var* must be a positive integer.
    :raises ValueError: Parameter *var* must be a positive integer
      smaller than nmax.
    '''
    if not isinstance(var, int):
        raise TypeError('\n\nParameter {} must be an integer.\n'.format(var_name))
    if var < 1:
        raise ValueError('\n\nParameter {} must be a positive integer.\n'.format(var_name))
    if var > nmax:
        raise ValueError('\n\nParameter {} must be a positive integer\n'\
                                   'smaller than {}.\n'.format(var_name, nmax))


def real_number(var, var_name):
    '''
    Check if parameter *var* is a real number.

    :raises TypeError: Parameter *var* must be a real number.
    '''         
    if not isinstance(var, (int, float)):
        raise TypeError('\n\nParameter {} must be a real number.\n'.format(var_name))


def positive_real(var, var_name):
    '''
    Check if parameter *var* is a positive number.

    :raises TypeError: Parameter *var* must be a real number.
    :raises ValueError: Parameter *var* must be a positive number.
    '''         
    if not isinstance(var, (int, float)):
        raise TypeError('\n\nParameter {} must be a real number.\n'.format(var_name))
    if var <= 0:
        raise ValueError('\n\nParameter {} must be a positive number.\n'.format(var_name))


def positive_real_zero(var, var_name):
    '''
    Check if parameter *var* is a positive number or zero.

    :raises TypeError: Parameter *var* must be a real number.
    :raises ValueError: Parameter *var* must be a positive number or zero.
    '''         
    if not isinstance(var, (int, float)):
        raise TypeError('\n\nParameter {} must be a real number.\n'.format(var_name))
    if var < 0:
        raise ValueError('\n\nParameter {} must be a positive number  or zero.\n'.format(var_name))


def negative_real(var, var_name):
    '''
    Check if parameter *var* is a negative number.

    :raises TypeError: Parameter *var* must be a real number.
    :raises ValueError: Parameter *var* must be a negative number.

    '''         
    if not isinstance(var, (int, float)):
        raise TypeError('\n\nParameter {} must be a real number.\n'.format(var_name))
    if var >= 0:
        raise ValueError('\n\nParameter {} must be a negative number.\n'.format(var_name))


def number(var, var_name):
    '''
    Check if parameter *var* is a number.

    :raises TypeError: Parameter *var* must be a real number.
    '''         
    if not isinstance(var, (int, float, complex)):
        raise TypeError('\n\nParameter {} must be a real number.\n'.format(var_name))


def larger(var1, var_name1, var2, var_name2):
    '''
    Check if *var1* larger than *val*.

    :raises ValueError: Parameter *var1* larger than *var2*.
    '''         
    if var1 >= var2:
        raise ValueError('\n\n{} must be larger than {}.\n'
                                    .format(var_name1, var_name2))


def smaller(var1, var_name1, var2, var_name2):
    '''
    Check if *var1* smaller than *var2*.

    :raises ValueError: Parameter *var1* must be smaller than *var2*.
    '''         
    if var1 >= var2:
        raise ValueError('\n\n{} must be smaller than {}.\n'
                                    .format(var_name1, var_name2))


def string(var, var_name):
    '''
    Check if parameter *var* is a string.

    :raises TypeError: Parameter *var* must be a string.
    '''
    if var is None:
        return
    if not isinstance(var, str):
        raise TypeError('\n\nParameter {} must be a string.\n'.format(var_name))


def ndarray(var, var_name, length):
    '''
    Check if parameter *var* is a numpy array.

    :raises TypeError: Parameter *var* must be a numpy ndarray.
    :raises ValueError: length array must be equal to length.
    '''
    if not isinstance(var, np.ndarray):
        raise TypeError('\n\nParameter {} must be a numpy ndarray.\n'.format(var_name))
    if len(var) != length:
        raise ValueError('\n\nParameter {} must be a numpy ndarray of {}\n'
                                    ''.format(var_name, length))


def ndarray_null(var, var_name):
    '''
    Check if parameter *var* is not a null numpy array.

    :raises ValueError: Parameter *var* must not be a null numpy ndarray.
    '''
    array_null = np.zeros(len(var))
    if np.allclose(var, array_null):
        raise ValueError('\n\nParameter {} must not be a null numpy ndarray.\n'.format(var_name))


def ndarray_empty(var, var_name):
    '''
    Check if parameter *var* is not an empty numpy array.

    :raises ValueError: Parameter *var* must not be an emptynumpy ndarray.
    '''
    if var.size == 0:
        raise ValueError('\n\nParameter {} must not be an emptynumpy ndarray.\n'.format(var_name))


def list_tuple_2elem(var, var_name):
    '''
    Check if parameter *var* is a list/tuple with 2 elements.

    :raises TypeError: Parameter *var* must be a list/tuple.
    :raises ValueError: Parameter *var* must contain 2 elements.
    '''
    if var is None:
        return
    if not isinstance(var, (list, tuple)):
        raise TypeError('\n\nParameter {} must be a list/tuple\n'.format(var_name))
    if len(var) != 2:
        raise ValueError('\n\nParameter {} must be a list/tuple of length two.\n'.format(var_name))


def tuple_2elem(var, var_name):
    '''
    Check if parameter *var* is a tuple with 2 elements.

    :raises TypeError: Parameter *var* must be a tuple.
    :raises ValueError: Parameter *var* must contain 2 elements.
    '''
    if var is None:
        return
    if not isinstance(var, tuple):
        raise TypeError('\n\nParameter {} must be a tuple\n'.format(var_name))
    if len(var) != 2:
        raise ValueError('\n\nParameter {} must be a list/tuple of length two.\n'.format(var_name))


###############################
# LATTICE EXCEPTION HANDLING
###############################


def lat(lat):
    '''
    Check if parameter is an instance of the *lattice*.
    :raises TypeError: Parameter must be an instance of the class lattice.
    '''
    names = inspect.getmro(lat.__class__)
    if str(names).find('lattice') == -1:
        raise TypeError('\n\nParameter must be an instance of the class lattice.\n')


def unit_cell(unit_cell):
    '''
    Check parameter *unit_cell*.

    :raises TypeError: Parameter unit_cell must be a list.
    :raises KeyError: Dictionaries must contain the key "tag".
    :raises KeyError: Dictionaries must contain the key "r0".
    :raises TypeError: Key "tags" must contain a binary char.
    :raises ValueError: Key "tags" must contain a binary char.
    :raises ValueError: Key "r0" must contain be a list.
    :raises TypeError: Key "r0" must contain be a tuple.
    :raises ValueError: Key "r0" must contain a tuple of length two.
    :raises ValueError: Key "r0" must contain a tuple of two real numbers.
    '''
    if not isinstance(unit_cell, list):
        raise TypeError('\n\nParameter unit_cell must be a list.\n')
    for dic in unit_cell:
        if 'tag' not in dic:
            raise KeyError('\n\nDictionaries must contain the key "tag".\n')
        if 'r0' not in dic:
            raise KeyError('\n\nDictionaries must contain the key "r0".\n')
        if not isinstance(dic['tag'], bytes):
            raise TypeError('\n\Key "tags" must contain a binary char.\n')    
        if not len(dic['tag']) == 1:
            raise ValueError('\n\Key "tags" must be a binary char.\n')    
        if not isinstance(dic['r0'], tuple):
            raise TypeError('\n\Key "r0" must be a tuple.\n')
        if not len(dic['r0']) == 2:
            raise ValueError('\n\Key "r0" must contain a tuple of length two.\n')    
        if not isinstance(dic['r0'][0], (int, float)) or not isinstance(dic['r0'][1], (int, float)):
            raise ValueError('\n\Key "r0" must contain a tuple of two real numbers.\n')

    
def prim_vec(prim_vec):
    '''
    Check parameter *prim_vec*.

    :raises TypeError: Parameter prim_vec must be a list.
    :raises ValueError: Parameter prim_vec must be a list
      of length 1 for 1D lattices or length 2 fro 2D lattices.
    :raises TypeError: List elements must be tuples.
    :raises ValueError: List elements must be 1 or 2 tuples.
    :raises ValueError: Tuples must be of length 2.
    :raises ValueError: Tuples must containt real numbers.
    :raises ValueError: Norm of coor should be larger than 0.1.
    '''
    if not isinstance(prim_vec, list):
        raise TypeError('\n\nParameter prim_vec must be a list.\n')
    if not len(prim_vec) == 1 and not len(prim_vec) == 2:
        raise ValueError('\n\nParameter prim_vec must be a list.\n'
                                  'of length 1 for 1D lattices or length 2 fro 2D lattices.\n')
    if 1 > len(prim_vec) > 2:
        raise ValueError('\n\nParameter prim_vec value must be a list of length 1 or 2.\n')
    for coor in prim_vec:
        if not isinstance(coor, tuple):
            raise TypeError('\n\nParameter prim_vec contain tuples\n')
        if len(coor) != 2:
            raise ValueError('\n\nParameter prim_vec contain tuples of length 2.\n')
        if (not isinstance(coor[0], (int, float))) or \
           (not isinstance(coor[1], (int, float))):
            raise ValueError('\n\ncoor must contain real numbers.\n')
        if coor[0] ** 2 + coor[1] ** 2 < 0.1:
            raise ValueError('\n\nNorm of coor should be larger than 0.1.\n')


def get_lattice(prim_vec, n1, n2):
    '''
    Check method *get_lattice*.

    :raises TypeError: Parameter n1 must be an integer.
    :raises TypeError: Parameter n2 must be an integer.
    :raises ValueError: Parameter n1 must be a positive integer.
    :raises ValueError: Parameter n2 must be a positive integer.
    '''
    positive_int(n1, 'n1')
    positive_int(n2, 'n2')
    if len(prim_vec) == 1 and n2 > 1:
        raise ValueError('\n\nParameter n2 should be equal to 1\n'
                                    'if one primitive vector is given\n')


def coor(coor):
    '''
    Check if *coor* is a structured array with 
    dtype=[('x', '<f8'), ('y', '<f8'), ('tag', 'S1')].
    '''
    if coor.dtype != [('x', '<f8'), ('y', '<f8'), ('tag', 'S1')]:
        raise TypeError('\n\nParameter coor dtype must be\n\
                                  dtype=[("x", "f8"), ("y", "f8"), ("tag", "S1")].\n')


def empty_coor(coor):
    '''
    Check if *get_lattice* has been called (*coor* not empty).
    :raises RuntimeError: Run method get_lattice first.
    '''
    if coor.size == 0:
        raise RuntimeError('\n\nRun method get_lattice first.\n')


def empty_coor_hop(coor):
    '''
    Check if *get_lattice_hop* has been called (*coor_hop* not empty).
    :raises RuntimeError: Run method sys.get_coor_hop first.
    '''
    if coor.size == 0:
        raise RuntimeError('\n\nRun method sys.get_coor_hop first.\n')


def coor_1d(coor):
    '''
    Check if *coor* is 1d (coor['y'] = cst).
    :raises ValueError: *coor* must be 1d( coor['y'] = cst)..
    '''
    if not np.allclose(coor['y'][0], coor['y']):
        raise ValueError('\n\ncoor["y"] must be constant.\n')


def remove_sites(index, sites):
    '''
    Check method *remove_sites*.

    :raises TypeError: Parameter index must be a list.
    :raises ValueError: Parameter index must be a list of integers.
    :raises ValueError: Indices must be between 0 and sites -1.
      of integers between 0 and sites
    '''
    if not isinstance(index, list):
        raise TypeError('\n\nParameter index must be a list.\n')
    if not all(isinstance(i, int) for i in index):
        raise ValueError('\n\nParameter index must be a list of integers.\n')
    if not all(-1 < i < sites for i in index):
        raise ValueError('\n\nElements of index must be between 0 and sites - 1.\n')


def shift(shift):
    '''
    Check *shift_x* and *shift_y*.
    :raises TypeError: Parameter delta must be a real number.
    '''
    if not isinstance(shift, (int, float)):
        raise TypeError('\n\nParameter shift must be a real number.\n')


def boundary_line(cx, cy, co):
    '''
    Check *boundary_line*.
    :raises TypeError: Parameter cx must be a real number.
    :raises TypeError: Parameter cy must be a real number.
    :raises TypeError: Parameter co must be a real number.
    '''
    if not isinstance(cx, (int, float)):
        raise TypeError('\n\nParameter cx must be a real number.\n')
    if not isinstance(cy, (int, float)):
        raise TypeError('\n\nParameter cy must be a real number.\n')
    if not isinstance(co, (int, float)):
        raise TypeError('\n\nParameter co must be a real number.\n')


def ellipse(a, b):
    '''
    Check *ellipse_in* and *ellipse_out*.
    :raises TypeError: Parameter a must be a positive number.
    :raises TypeError: Parameter b must be a positive number.
    '''
    if not isinstance(a, (int, float)):
        raise TypeError('\n\nParameter a must be a positive number.\n')
    if not isinstance(b, (int, float)):
        raise TypeError('\n\nParameter b must be a positive number.\n')
    if a <= 0:
        raise ValueError('\n\nParameter a must be a positive number.\n')
    if b <= 0:
        raise ValueError('\n\nParameter b must be a positive number.\n')


def sites(sites):
    '''
    Check if *get_lattice* has been called (*coor* not empty).
    :raises RuntimeError: Run method lat.get_lattice first.
    '''
    if sites == 0:
        raise RuntimeError('\n\nRun method lat.get_lattice first.\n')


####################################
# CLASS SYSTEM EXCEPTION HANDLING
####################################


def sys(sys):
    '''
    Check if parameter is an instance of the *system*.
    :raises TypeError: Parameter must be an instance of the class system.
    '''
    names = inspect.getmro(sys.__class__)
    if str(names).find('system') == -1:
        raise TypeError('\n\nParameter must be an instance of the class system.\n')


def print_hopping(n, nmax):
    '''
    Check method *print_vec_hopping*.

    :raises TypeError: Parameter *nmax* must be an integer.
    :raises ValueError: Parameter *nmax* must be a positive integer.
      between 1 and n_max-1.
    '''
    if not isinstance(n, int):
        raise TypeError('\n\nParameter n_max must be an integer.\n')
    if n < 1 or n > nmax-1:
        raise ValueError('\n\nParameter n_max must be a positive integer'
                                    'between 1 and n_max-1.\n')


def set_onsite(onsite, tags):
    '''
    Check method *set_onsite*.
    :raises TypeError: Parameter onsite must be a dictionary.
    :raises ValueError: Parameter onsite keys must be a tag.
    :raises ValueError: Parameter onsite values must be
      real and/or complex numbers.
    '''
    if not isinstance(onsite, dict):
        raise TypeError('\n\nParameter onsite must be a dictionary.\n')
    for tag, val in onsite.items():
        if tag not in tags:
            raise ValueError('\n\nParameter onsite keys must be a tag.\n')   
        if not isinstance(val, (int, float, complex)):
            raise ValueError('\n\nParameter onsite values must be\n'\
                                       'real and/or complex numbers.\n')


def set_hopping(list_hop, n_max):
    '''
    Check method *set_hopping*.

    :raises TypeError: Parameter *list_hop* must be a list.
    :raises TypeError: Parameter *list_hop* must be a list of dictionary.
    :raises KeyError: "n" and "t" must be dictionary keys.
    :raises KeyError: "tag" or "ang" must be a key.
    :raises KeyError: "tag" and "ang" must be a key.
    :raises ValueError: Dictionaries must be of length 2, 4, or 4.
    :raises ValueError: "n" must be between 1 and nmax"

    '''
    if not isinstance(list_hop, list):
        raise TypeError('\n\nParameter *list_hop* must be a list.\n')
    for dic in list_hop:
        if not isinstance(dic, dict):
            raise TypeError('\n\nParameter *list_hop* must be a list of dictionary.\n')
        if 'n' not in dic or 't' not in dic:
                raise KeyError('\n\n"n" and "t" must be dictionary keys.\n')
        if not isinstance(dic['n'], int):
            raise TypeError('\n\n"n" value must be an integer.\n')
        if not 0 < dic['n'] <= n_max:
            raise ValueError('\n\n"n" value must be between 1 and nmax".\n')
        if not isinstance(dic['t'], (int, float, complex)):
            raise TypeError('\n\n"t" value must be a real or complex number.\n')
        if len(dic) == 3:
            if not ('tag' not in dic or 'ang' not in dic):
                raise KeyError('\n\n"tag" or "ang" must be a key.\n')
        elif len(dic) == 4:
            if 'tag' not in dic and 'ang' not in dic:
                raise KeyError('\n\n"tag" or "ang" must be a key.\n')
        elif len(dic) > 4:
            raise ValueError('\n\nDictionaries must be of length 2, 3, or 4.\n')
        if 'tag' in dic:
            if not isinstance(dic['tag'], bytes):
                raise TypeError('\n\n"tag" value must be a binary string.\n')
            if len(dic['tag']) != 2:
                raise ValueError('\n\n"tag" value must be a binary string of length 2.\n')
        if 'ang' in dic:
            if not isinstance(dic['ang'], (int, float)):
                raise TypeError('\n\n"ang" value must be a real number.\n')


def index(ind, dic):
    '''
    check if *ind* not empy.
    '''
    if np.sum(ind) == 0:
        raise ValueError('\n\nNo hoppings with parameters {}.\n'.format(dic))


def set_hopping_def(hop, hopping_def, sites):
    '''
    Check method *set_hop_def*.

    :raises TypeError: Parameter *hopping_def* must be a dictionary
    :raises TypeError: *hopping_def* keys must be lists.
    :raises ValueError: *hopping_def* keys must be lists of length 2.
    :raises ValueError: *hopping_def* keys must be lists of integers.
    :raises TypeError: *hopping_def* keys must be lists.
    :raises ValueError: *hopping_def* keys must be integers between 0 and sites-1.
    :raises ValueError: *hopping_def* keys must be different integers between 0 and sites-1.
    :raises TypeError: *hopping_def* values must be numbers.
    '''
    if not isinstance(hopping_def, dict):
        raise TypeError('\n\nParameter hopping_def must be a dictionary.\n')
    for key, val in hopping_def.items():
        if not isinstance(key, tuple):
            raise TypeError('\n\nhopping_def keys must be lists.\n')
        if len(key) != 2:
            raise TypeError('\n\nhopping_def keys must be lists of length 2.\n')
        if not isinstance(key[0], int) or not isinstance(key[1], int):
            raise ValueError('\n\nhopping_def keys must be lists of integers.\n')
        if key[0] < 0 or key[1] < 0 or key[0] > sites-1 or key[1] > sites-1:
            raise ValueError('\n\nhopping_def keys must be integers between 0 and sites-1.\n')
        if key[0] == key[1]:
            raise ValueError('\n\nhopping_def keys must be different integers between 0 and sites-1.\n')
        if not isinstance(val, (int, float, complex)):
            raise TypeError('\n\nhopping_def values must be numbers.\n')


def set_onsite_def(onsite_def, sites):
    '''
    Check method *set_ons_def*.

    :raises TypeError: Parameter *onsite_def* must be a dictionary.
    :raises TypeError: *onsite_def* keys must be integers.
    :raises TypeError: *onsite_def* values must be numbers.
    :raises ValueError: *onsite_def* keys must be integers between :math:`[0, sites)`.
    '''
    if not isinstance(onsite_def, dict):
        raise TypeError('\n\nParameter onsite_def must be a dictionary.\n')
    for key, val in onsite_def.items():
        if not isinstance(key, int):
            raise TypeError('\n\nonsite_def keys must be integers.\n')
        if not isinstance(val, (int, float, complex)):
            raise TypeError('\n\nonsite_def values must be numbers.\n')
        if key < 0 or key > sites-1:
            raise ValueError('\n\nonsite_def keys must be integers between 0 and sites-1.\n')


def hop_n1(hop):
    '''
    Check method if self.hop contains nearest neighbours hoppings.

    :raises RunTimeError: self.hop must contain nearest neighbours hoppings.
    '''
    if len(hop['n'] == 1) == 0:
        raise ValueError('\n\nParameter hop must contain nearest neighbours hoppings.\n')


def empty_onsite(onsite):
    '''
    Check if *onsite* not empty.

    :raises RuntimeError: Run method set_onsite first.
    '''
    if onsite.size == 0:
        raise RuntimeError('\n\nRun method set_onsite first\n')


def empty_hop(hop):
    '''
    Check if *hop* not empty.

    :raises RuntimeError: Run method set_hopping first.
    '''
    if hop.size == 0:
        raise RuntimeError('\n\nRun method set_hopping first\n')


def hop_sites(hop, sites):
    '''
    Check if *hop* indices are smaller than *sites*.

    :raises ValueError: Run method sys.clean_hopping.
    '''
    row_max = np.max(hop['i'])
    col_max = np.max(hop['i'])
    ind_max = max(row_max, col_max)
    if sites < ind_max:
        raise ValueError('\n\nRun method system.clean_hopping.\n'
                                    'and redefine the hoppings.\n')


def empty_coor(coor):
    '''
    Check if *coor* not empty.

    :raises RuntimeError: Run method lattice.get_lattice first.
    '''
    if coor.size == 0:
        raise RuntimeError('\n\nRun method lattice.get_lattice first.\n')


def empty_coor_hop(coor_hop):
    '''
    Check if *coor_hop* not empty.

    :raises RuntimeError: Run method system.get_coor_hop first.
    '''
    if coor_hop.size == 0:
        raise RuntimeError('\n\nRun method system.get_coor_hop first.\n')


def empty_ham(ham):
    '''
    Check if Hamiltonian not empty.

    :raises RuntimeError: Run method system.get_ham first.
    '''
    if not ham.nnz:
        raise RuntimeError('\n\nRun method system.get_ham first.\n')


def empty_en(en):
    '''
    Check if *en* not empty.

    :raises RuntimeError: Run method get_ham first.
    '''
    if en.size == 0:
        raise RuntimeError('\n\nRun method get_eig first\n')


def empty_pola(pola):
    '''
    Check if *pola* not empty.

    :raises RuntimeError: Run method get_eig(eigenvec=True) first.
    '''
    if pola.size == 0:
        raise RuntimeError('\n\nRun method get_eig(eigenvec=True) first\n')


def empty_vn(vn):
    '''
    Check if *vn* not empty.

    :raises RuntimeError: Run method get_eig(eigenvec=True) first.
    '''
    if vn.size == 0:
        raise RuntimeError('\n\nRun method get_eig(eigenvec=True) first\n')


def empty_ipr(ipr):
    '''
    Check if *ipr* not empty.
    '''
    if ipr.size == 0:
        raise RuntimeError('\n\nRun method get_ipr first\n')


def empty_ndarray(arr, method):
    '''
    Check if *arr* is a not empty np.ndarray.
    '''
    if arr.size == 0:
        raise RuntimeError('\n\nRun method {} first\n'. format(method))


def tag(tag, tags):
    '''
    Check tag.

    :raises TypeError: Parameter *tag* must be a binary string.
    :raises ValueError: Parameter *tag* is not in tags.
    '''
    if not isinstance(tag, bytes):
        raise TypeError('\n\nParameter tag must be a binary char.\n')
    if tag not in tags:
        raise ValueError('\n\nParameter tag is not in tags.\n')


def angle(angle, angles, upper_part):
    '''
    Check angle.

    :raises TypeError: Parameter *angle* must be

        * a positive number if *upper_part* is True
        * a negative real if *upper_part* is False.
    :raises ValueError: Parameter *angle* is not in hop['ang'].
    :raises ValueError: Parameter *angle*must be positive.
    '''
    if upper_part:
        positive_real_zero(angle, 'angle, if upper_part=True,')
    else:
        negative_real(angle, 'angle, if upper_part=False,')
    if not np.sum(np.isclose(angle, angles)) and np.sum(np.isclose(angle, 180+angles)):
        raise ValueError('\n\nParameter angle is not in hop["ang"].\n')
    

def lims(lims):
    '''
    Check parameter *lims*.

    :raises TypeError: Parameter lims must be a list.
    :raises TypeError: Parameter *lims[0]* must be a real number.
    :raises TypeError: Parameter *lims[1]* must be a real number.
    :raises ValueError: *lims* must be a list of length 2.
    :raises ValueError: *lims[0]* must be smaller than *lims[1]*.
    '''
    if lims is not None:
        list_tuple_2elem(lims, 'lims')
        real_number(lims[0], 'lims[0]')
        real_number(lims[1], 'lims[1]')
        smaller(lims[0], 'lims[0]', lims[1], 'lims[1]')


def lims_positive(lims):
    '''
    Check parameter *lims*.

    :raises TypeError: Parameter lims must be a list.
    :raises TypeError: Parameter *lims[0]* must be a positive real number.
    :raises TypeError: Parameter *lims[1]* must be a positive real number.
    :raises ValueError: *lims* must be a list of length 2.
    :raises ValueError: *lims[0]* must be smaller than *lims[1]*.
    '''
    if lims is not None:
        list_tuple_2elem(lims, 'lims')
        positive_real(lims[0], 'lims[0]')
        positive_real(lims[1], 'lims[1]')
        smaller(lims[0], 'lims[0]', lims[1], 'lims[1]')


#################################
# CLASS PLOT EXCEPTION HANDLING
#################################


def fig(fig):
    '''
    Check if fig is an instance of *Figure*.

    :raises TypeError: fig must be an instance of *Figure*.
    '''
    if not fig.__class__.__name__ == 'Figure':
        raise TypeError('\n\nfig must be an instance of *Figure*.\n')


def ani(ani):
    '''
    Check if ani is an instance of *FuncAnimation*.

    :raises TypeError: ani must be an instance of *FuncAnimation*.
    '''
    if not fig.__class__.__name__ == 'FuncAnimation':
        raise TypeError('\n\nani must be an instance of *FuncAnimation*.\n')


def file_format(file_format):
    '''
    Check if file_format is a string 'png', 'pdf', 'ps', 'eps', or 'svg'.

    :raises TypeError: file_format must be a string.
    :raises ValueError: file_format must be a string given by 'png', 'pdf', 'ps', 'eps', or 'svg'.

    '''
    if not isinstance(file_format, str):
        raise TypeError('\n\nfile_format must be a string.\n')
    if file_format not in ['png', 'pdf', 'ps', 'eps', 'svg']:
        raise ValueError('\n\nfile_format must be a string given by,\n'\
                                   ' "png", "pdf", "ps", "eps", or "svg".\n')


def fig_size(sigsize):
    if figsize is None:
        figsize = (5, 4)
    error_handling.list_tuple_2elem(figsize, 'figsize')
    error_handling.positive_real(figsize[0], 'figsize[0]')
    error_handling.positive_real(figsize[1], 'figsize[1]')


####################################
# PROPAGATION
####################################


def get_pump(hams):
    if not isinstance(hams, list):
        raise TypeError('\n\nhams must be a list.\n')
    for ham in hams:
        empty_ham(ham)

def prop_type(prop_type):
    string(prop_type, 'prop_type')
    if prop_type not in ['real', 'imag', 'norm']:
        raise ValueError('\n\nParameter prop_type must be a string:\n'
                                   '"real", "imag", "norm".\n')
