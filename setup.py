"""Set up MEICA module.
"""
from setuptools import setup

# fetch version from within meica module
with open('meica/version.py') as f:
    exec(f.read())

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
    entry_points={'console_scripts': [
        'meica=meica.meica:main'
    ]},
    license='LGPLv2+',
    **extra_setuptools_args)
