from setuptools import setup, find_packages


packages=find_packages(exclude=['contrib', 'docs', 'tests']),
print(packages)
setup(  name = 'tbee', 
            version = '0.1',
            packages=find_packages(exclude=['contrib', 'docs', 'tests']),
            test_suite = 'tests',
        )
DESCRIPTION = "Tight-Binding module"
LONG_DESCRIPTION = "Tight-Binding module"
NAME = "tbee"
AUTHOR = "Charles Poli"
AUTHOR_EMAIL = "cpoli83@hotmail.fr"
MAINTAINER = "Charles Poli"
MAINTAINER_EMAIL = "jcpoli83@hotmail.fr"
DOWNLOAD_URL = 'https://github.com/cpoli/tbee'
LICENSE = 'BSD'

