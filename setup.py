from setuptools import setup, find_packages


packages=find_packages(exclude=['contrib', 'docs', 'tests']),
print(packages)
setup(  name = 'TB', 
            version = '0.1',
            #packages = find_packages(),
            packages=find_packages(exclude=['contrib', 'docs', 'tests']),
            test_suite = 'tests',
            #install_requires=[]
        )
DESCRIPTION = "Tight-Binding module"
LONG_DESCRIPTION = DESCRIPTION
NAME = "TB"
AUTHOR = "Charles Poli"
AUTHOR_EMAIL = "cpoli83@hotmail.fr"
MAINTAINER = "Charles Poli"
MAINTAINER_EMAIL = "jcpoli83@hotmail.fr"
DOWNLOAD_URL = 'https://github.com/cpoli/TB'
LICENSE = 'BSD'

print(packages)
'''
setup(name=NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      maintainer=MAINTAINER,
      maintainer_email=MAINTAINER_EMAIL,
      url=DOWNLOAD_URL,
      download_url=DOWNLOAD_URL,
      license=LICENSE,
      packages=['TB'],
#      package_data={'JSAnimation': ['icons/*.png']}
     )
'''