from __future__ import annotations

import os
from matplotlib.figure import Figure
from matplotlib.animation import FuncAnimation
import tbee.error_handling as error_handling


#################################
# CLASS SAVE
#################################


class Save():
    '''
    Create folder and save figures / animations obtained via
    **plot** or **propagation**.

    :param dir_name: String. Name of the sub-directory the figures are stored in.
    :param dir_main: String. Default value None. Path of the main directory.
      If None, figures are stored under ``'figs/'``.
    :param params: Dictionary. Default value None. Parameters appended to file names.
    :param file_format: Default value 'png'. Figure format.
    '''
    def __init__(
        self,
        dir_name: str,
        dir_main: str | None = None,
        params: dict | None = None,
        file_format: str = 'png',
    ) -> None:
        error_handling.string(dir_name, 'dir_name')
        error_handling.string(dir_main, 'dir_main')
        error_handling.file_format(file_format)
        self.params = {} if params is None else params
        self.file_format = file_format
        if dir_main is None:
            self.dir_main = 'figs/'
        else:
            self.dir_main = dir_main
        self.dir_name = self.dir_main + dir_name
        self.create_dir()

    def create_dir(self) -> None:
        '''
        Create the directory to store the figures exists.
        '''
        if not os.path.exists(self.dir_main):
            os.makedirs(self.dir_main)
        if not os.path.exists(self.dir_name):
            os.makedirs(self.dir_name)

    def file_name(self) -> str:
        '''
        Create the file name.

        :returns:
            * **file_name** -- File name.
        '''
        file_name = ''
        for key, val in self.params.items():
            file_name += '_' + key + str(complex(val+0)).replace('.', ',')
        return file_name

    def fig(self, fig: Figure, name: str) -> None:
        '''
        Save the figure in the directory defined by the method *dir_name()*.

        :param fig: Matplotlib fig.
        :param name:  String. Fist part of the file name.
        '''
        error_handling.fig(fig)
        error_handling.string(name, 'name')
        name_file = self.dir_name + '/' + name + self.file_name() + '.' + self.file_format
        fig.savefig(name_file, format=self.file_format)

    def fig_lat(self, fig: Figure, name: str) -> None:
        '''
        Save the figure in the directory defined by the method *dir_name()*.

        :param fig: Matplotlib fig.
        :param name:  String. First part of the file name.
        '''
        error_handling.fig(fig)
        error_handling.string(name, 'name')
        name_file = self.dir_name + '/' + name + '.' + self.file_format
        fig.savefig(name_file, format=self.file_format)

    def ani(self, ani: FuncAnimation, name: str, fps: int = 10) -> None:
        error_handling.ani(ani)
        error_handling.string(name, 'name')
        error_handling.positive_int(fps, 'fps')
        name_file = self.dir_name + '/' + name + '.mp4'
        ani.save(name_file, fps=fps, extra_args=['-vcodec', 'libx264'])


# Backward-compatible lowercase alias (pre-0.2 API).
save = Save
