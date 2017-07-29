from setuptools import setup

__version__ = '3.2.1a'

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

extra_setuptools_args = dict(
    tests_require=['pytest'])

setup(
    name='meqc',
    version=__version__,
    description='Multi-echo independent components analysis of fMRI data.',
    maintainer='Prantik Kundu',
    maintainer_email='',
    install_requires=requirements,
    packages=['meqc'],
    license='MIT',
    **extra_setuptools_args)
