from setuptools import setup, find_packages

setup(
    name='rodriguez-lab-notebook-franklinhiciano',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
      'biopython',
      'pandas',
      'receptor-utils',
    ],
    author='Franklin Hiciano',
    author_email='fhiciano5@gmail.com',
    description="Franklin Hiciano's notebook for Rodriguez Lab at Mount Sinai",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/franklin-hiciano/rodriguez-lab-notebook',
    classifiers=[
      'Programming Language :: Python :: 3',
      'License :: OSI Approved :: MIT License',
      'Operating System :: OS Independent',
    ],
)
