from tbee.save import Save, save
import unittest
import os
import shutil
import tempfile
import matplotlib.pyplot as plt


class TestSave(unittest.TestCase):
    '''
    Unittest of class **Save**.
    '''
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.cwd = os.getcwd()
        os.chdir(self.tmp)

    def tearDown(self):
        os.chdir(self.cwd)
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_default_dir_main(self):
        sav = Save(dir_name='test1')
        self.assertEqual(sav.dir_main, 'figs/')
        self.assertTrue(os.path.exists(sav.dir_name))

    def test_explicit_dir_main(self):
        sav = Save(dir_name='test2', dir_main='other/')
        self.assertEqual(sav.dir_main, 'other/')
        self.assertTrue(os.path.exists(sav.dir_name))

    def test_create_dir_idempotent(self):
        sav = Save(dir_name='test3')
        sav.create_dir()  # directories already exist: no-op branch.
        self.assertTrue(os.path.exists(sav.dir_name))

    def test_file_name(self):
        sav = Save(dir_name='test4', params={'t': 1.5, 'e': -1j})
        name = sav.file_name()
        self.assertTrue(name.startswith('_'))

    def test_fig(self):
        sav = Save(dir_name='test5', params={'t': 1.})
        fig = plt.figure()
        sav.fig(fig, 'spectrum')
        plt.close(fig)
        files = os.listdir(sav.dir_name)
        self.assertEqual(len(files), 1)
        self.assertTrue(files[0].startswith('spectrum'))
        self.assertRaises(TypeError, sav.fig, 0, 'spectrum')

    def test_fig_lat(self):
        sav = Save(dir_name='test6')
        fig = plt.figure()
        sav.fig_lat(fig, 'lattice')
        plt.close(fig)
        self.assertTrue(os.path.exists(os.path.join(sav.dir_name, 'lattice.png')))
        self.assertRaises(TypeError, sav.fig_lat, 0, 'lattice')

    def test_ani(self):
        class FuncAnimation:
            def save(self, name_file, fps, extra_args):
                with open(name_file, 'w') as f:
                    f.write('fake animation')

        sav = Save(dir_name='test7')
        sav.ani(FuncAnimation(), 'prop', fps=5)
        self.assertTrue(os.path.exists(os.path.join(sav.dir_name, 'prop.mp4')))
        self.assertRaises(TypeError, sav.ani, 0, 'prop')

    def test_backward_compatible_alias(self):
        self.assertIs(save, Save)


if __name__ == '__main__':
    unittest.main()
