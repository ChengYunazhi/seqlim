from distutils.core import setup
from seqlim import __version__

setup(
    name='seqlim',
    version=__version__,
    packages=['seqlim',],
    license='Apache 2.0',
    long_description=open('README.md').read(),
    scripts = ['bin/seqlim']
)
