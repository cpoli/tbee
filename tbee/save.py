import os
import tbee.error_handling as error_handling


#################################
# CLASS SAVE
#################################


class save():
    '''
    Create folder and save figures / animations obtained via 
    **plot** or **propagation**. 
    Plot the results of the classes **lattice** or **system**.

    :param dir_main: Name of the directory.
    :param dir_name: Default value None. Relative path of the main directory.
      if None, figures stored in ''../TBfig/'dir_name/'
    :param params: dictionary. file name information
    :param file_format: Default value 'png'. Figure format.
    '''
    def __init__(self, dir_name, dir_main=None, params={}, file_format='png'):
        error_handling.string(dir_name, 'dir_name')
        error_handling.string(dir_main, 'dir_main')
        error_handling.file_format(file_format)
        self.params = params
        self.file_format = file_format
        if dir_main is None:
            self.dir_main = 'figs/'
        else:
            self.dir_main = dir_name
        self.dir_name = self.dir_name(dir_name)
        self.create_dir()

    def dir_name(self, dir_name):
        '''
        Set the name of the directory in which the figures are stored.

        :param dir_name: String. Directory name. 
        '''
        error_handling.string(dir_name, 'dir_name')
        dir_name = self.dir_main + dir_name
        return dir_name

    def create_dir(self):
        '''
        Create the directory to store the figures exists.
        '''
        if not os.path.exists(self.dir_main):
            os.makedirs(self.dir_main)
        if not os.path.exists(self.dir_name):
            os.makedirs(self.dir_name)

    def file_name(self):
        '''
        Create the file name.

        :returns:
            * **file_name** -- File name.
        '''
        file_name = ''
        for key, val in self.params.items():
            file_name += '_' + key + str(complex(val+0)).replace('.', ',')
        return file_name

    def fig(self, fig, name):
        '''
        Save the figure in the directory defined by the method *dir_name()*.

        :param fig: Matplotlib fig.
        :param name:  String. Fist part of the file name.
        '''
        error_handling.fig(fig)
        error_handling.string(name, 'name')
        name_file = self.dir_name + '/' + name + self.file_name() + '.' + self.file_format
        fig.savefig(name_file, format=self.file_format)

    def fig_lat(self, fig, name):
        '''
        Save the figure in the directory defined by the method *dir_name()*.

        :param fig: Matplotlib fig.
        :param name:  String. First part of the file name.
        '''
        error_handling.fig(fig)
        error_handling.string(name, 'name')
        name_file = self.dir_name + '/' + name + '.' + self.file_format
        fig.savefig(name_file, format=self.file_format)

    def ani(self, ani, name, fps=10):
        error_handling.ani(ani)
        error_handling.string(name, 'name')
        error_handling.positive_int(fps, 'fps')
        name_file = self.dir_name + '/' + name + '.mp4'
        ani.save(name_file, fps=fps, extra_args=['-vcodec', 'libx264'])
