#!/Users/suwang/miniconda3/bin/python
#!/usr/bin/python

import os
import sys

from setuptools import setup
from pkg_resources import resource_filename
from subprocess import call as subpcall
from setuptools import find_packages

if sys.version < "2.6.0":
    print("Please use a Python with higher version than 2.6.0")
    sys.stderr.write("CRITICAL: Python version must be higher than 2.6.0!\n")
    sys.exit(1)
    exit(1)
    
    
def main():

    #compilemis()
    if not float(sys.version[:3])>=2.6:
        sys.stderr.write("CRITICAL: Python version must be greater than or equal to 2.6! python 2.6 or newer is recommended!\n")
        sys.exit(1)
    setup(name="marge",
          version="1.0",
          description="MARGE -- Model-based Analysis of Regulation of Gene Expression ",
          author='Su Wang',
          author_email='wangsu0623@gmail.com',
          install_requires=['argparse','numpy','six','twobitreader','sklearn', 'jinja2'],
          packages=['marge'],
          package_data = {'marge': ['Snakefile', 'config.json']},
          zip_safe=False,
          scripts=['marge/marge'],
          classifiers=[
            'Development Status :: 4 - Beta',
            'Environment :: Console',
            'Environment :: Web Environment',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: Artistic License',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: POSIX',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            ],
          )
    
    
if __name__ == '__main__':
    main()
    
