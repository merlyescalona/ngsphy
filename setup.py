from setuptools import setup, find_packages
def readme():
    with open('README.md') as f:
        return f.read()
long_description='NGSphy is a Python open-source tool for the genome-wide simulation of NGS data (read counts or Illumina reads) obtained from thousands of gene families evolving under a common species tree, with multiple haploid and/or diploid individuals per species, where sequencing coverage (depth) heterogeneity can vary among species, individuals and loci, including off-target or uncaptured loci.'

setup(name='ngsphy',
    version='1.0.1',
    description='phylogenomic simulation of NGS data ',
    long_description=long_description,
    url='https://github.com/merlyescalona/ngsphy',
    download_url='https://github.com/merlyescalona/ngsphy/blob/master/dist/ngsphy-1.0.0.tar.gz',
    author='Merly Escalona',
    author_email='merlyescalona@uvigo.es',
    license='GNU/GPL v3',
    packages=['ngsphy'],
    package_dir={'ngsphy': 'ngsphy'},
    py_modules = [\
        'ngsphy.coverage', \
        'ngsphy.individual', \
        'ngsphy.loggingformatter', \
        'ngsphy.msatools', \
        'ngsphy.readcounts', \
        'ngsphy.reads', \
        'ngsphy.rerooter', \
        'ngsphy.sequence', \
        'ngsphy.settings'
    ],
    install_requires=[
        'argparse',\
        'ConfigParser',\
        'datetime',\
        'dendropy',\
        'logging',\
        'multiprocessing',\
        'numpy',\
        'setuptools',\
        'scipy'
    ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 2.7'
    ],
    keywords='biology phylogenomics next-generation sequencing coverage targeted-sequencing',
    python_requires='~=2.7',
    scripts=['scripts/ngsphy'],\
    entry_points={
        'console_scripts':[\
            'ngsphy = ngsphy.__main__:main'\
        ]\
    },
    zip_safe=False
  )
