#!/usr/bin/env python
# -*- coding: utf-8 -*-

from wiz.version import __version__

from setuptools import setup, find_packages


url = 'https://github.com/HadrienG/wiz'

with open('README.md') as f:
    long_description = f.read()

setup(
    name='wiz',
    version=__version__,

    description='takes a metagenome-assembled-genome and does a lot of fancy stuff with it',
    long_description=long_description,
    long_description_content_type='text/markdown',

    url=url,
    download_url=url + '/tarball/' + __version__,
    author='Hadrien Gourlé',
    author_email='hadrien.gourle@slu.se',

    license='MIT',
    packages=find_packages(),

    tests_require=['pytest, pytest-cov'],
    install_requires=['requests', 'tqdm', 'biopython']
    include_package_data=True,

    entry_points={
        'console_scripts': ['wiz = wiz.app:main'],
    }
)
