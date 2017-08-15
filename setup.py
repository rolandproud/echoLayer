# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 19:31:23 2017

@author: Roland Proud
"""

from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()
    
setup(name='pyechomask',
      version='1.0.0.dev2',
      description='An echogram classification tool',
      long_description=readme(),
      classifiers=[
        'Development Status :: 3 - Alpha',
        #'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        #'Topic :: Text Processing :: Linguistic',
      ],
      keywords='echosounder marine acoustics mask',
      url='https://github.com/rolandproud/pyechomask',
      author='Roland Proud',
      author_email='rp43@st-andrews.ac.uk',
      license='GPLv3',
      packages=['pyechomask'],
            install_requires=[
          'numpy',
          'matplotlib'
      ],
      python_requires='>=3',
      include_package_data=True,
      zip_safe=False)