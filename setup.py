"""A setuptools based setup module.
See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='cafa',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version='1.0.8',

    description='CAFA evaluator code',
    
    long_description='Evaluation code used for the 5th iteration of the CAFA challenge',

    # The project's main homepage.
    url='https://github.com/BioComputingUP/CAFA-evaluator',

    # Author details
    author='Damiano Piovesan',
    author_email='damiano.piovesan@unipd.it, parnal@iastate.edu',

    # Choose your license
    license='GPLv3',
    
    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
       # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],

    # What does your project relate to?
    keywords='CAFA evaluation',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),
    include_package_data=True,

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=['pandas', 'matplotlib', 'numpy'],

    # package_dir={'': 'cafa'},
    entry_points={
        'console_scripts': [
            'cafa = cafa.main:main',    # console script can be changed 
        ]
    },
)