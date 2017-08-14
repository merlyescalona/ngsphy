from setuptools import setup, find_packages
def readme():
    with open('README.md') as f:
        return f.read()

setup(name='ngsphy',
    version='1.0.0',
    description='phylogenomic simulation of NGS data ',
    long_description='',
    url='https://github.com/merlyescalona/ngsphy',
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
        'datetime',\
        'ConfigParser',\
        'datetime',\
        'dendropy',\
        'logging',\
        'multiprocessing',\
        'numpy',\
        'scipy'
    ],
    classifiers=[
        'Develpment Status :: 3 - Alpha',
        'Programming Language :: Python :: 2.7'
    ],
    keywords='biology phylogenomics next-generation sequencing ',
    python_requires='~=2.7',
    scripts=['scripts/ngsphy'],\
    entry_points={
        'console_scripts':[\
            'ngsphy = ngsphy.__main__:main'\
        ]\
    },
    zip_safe=False
  )
