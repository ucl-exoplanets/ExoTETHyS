#!/usr/bin/env python

"""
Setup file
"""

import io
import os
import re

from setuptools import setup, find_packages


def read(*names, **kwargs):
    with io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8")
    ) as fp:
        return fp.read()


def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


setup(
    name="exotethys",
    version=find_version("exotethys", "__init__.py"),
    description="",
    url="https://github.com/ucl-exoplanets/ExoTETHyS",
    author='Giuseppe Morello',
    author_email='giuseppe.morello.11@ucl.ac.uk',
    license='GNU General Public License v3 (GPLv3)',
    classifiers=['Development Status :: 4 - Beta',
                 'Environment :: Console',
                 'Intended Audience :: Science/Research',
                 'Topic :: Scientific/Engineering :: Astronomy',
                 'Operating System :: MacOS :: MacOS X',
                 'Programming Language :: Python :: 3.7',
                 ],
    packages=find_packages(),
    install_requires = ['scipy>=1.2.2', 'numpy>=1.16.5'],
    package_data={"exotethys":["*.pickle", "Passbands/*.pass"]},
    scripts=[],
    data_files=[('', ['README.md'])]
)
