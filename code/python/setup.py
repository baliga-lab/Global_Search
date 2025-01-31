import os
import re
import codecs
from setuptools import setup

NAME = 'globalsearch'
PACKAGES = ['globalsearch.rnaseq', 'globalsearch.control']
DESCRIPTION = 'globalsearch is a collection of Python modules and command tools for the Global Search pipeline.'
LICENSE = 'LGPL V3'
URI = 'https://github.com/baliga-lab/Global_Search'
AUTHOR = 'Wei-ju Wu'
VERSION = '0.3.0'

KEYWORDS = ['global search', 'coral reef', 'rna sequencing']

# See trove classifiers
# https://testpypi.python.org/pypi?%3Aaction=list_classifiers

CLASSIFIERS = [
    "Development Status :: 4 - Beta",
    "Environment :: Console",
    "Intended Audience :: Developers",
    "Natural Language :: English",
    "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    "Operating System :: OS Independent",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: Implementation :: CPython",
    "Topic :: Software Development :: Libraries :: Python Modules"
    ]
INSTALL_REQUIRES = ['jinja2', 'fs', 'rpy2']


if __name__ == '__main__':
    setup(name=NAME, description=DESCRIPTION,
          license=LICENSE,
          url=URI,
          version=VERSION,
          author=AUTHOR,
          author_email='weiju.wu@gmail.com',
          maintainer=AUTHOR,
          maintainer_email='weiju.wu@gmail.com',
          keywords=KEYWORDS,
          packages=PACKAGES,
          zip_safe=False,
          classifiers=CLASSIFIERS,
          install_requires=INSTALL_REQUIRES,
          scripts=['bin/gs_submit', 'bin/gs_prepare', 'bin/gs_run'])
