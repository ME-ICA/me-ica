from setuptools import setup

__version__ = '3.2.1'

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

extra_setuptools_args = dict(
    tests_require=['pytest'])

setup(
    name='meica',
    version=__version__,
    description='Multi-echo independent components analysis of fMRI data.',
    maintainer='Elizabeth DuPre',
    maintainer_email='emd222@cornell.edu',
    install_requires=requirements,
    packages=['meica'],
    license='LGPLv2+',
    **extra_setuptools_args)
