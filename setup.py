from setuptools import setup, find_packages
from gget.__init__ import __version__, __author__, __email__

def read(path):
    with open(path, 'r') as f:
        return f.read()

long_description = read('README.md')

setup(
    name='gget',
    version=__version__,
    license='MIT',
    author=__author__,
    author_email=__email__,
    maintainer=__author__,
    maintainer_email=__email__,
    description='Query Ensembl for genes using free form search words, look up genes/transcripts by Ensembl ID or fetch the latest FTPs by species.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    zip_safe=False,
    packages=find_packages(),
    include_package_data=True,
    python_requires='>=3.6',
    install_requires=read('requirements.txt').strip().split('\n'),
    url='https://github.com/lauraluebbert/gget',
    keywords='gget',
    entry_points={
        'console_scripts': ['gget=gget.main:main'],
    },
    classifiers=[
        'Environment :: Console',
        'Framework :: Jupyter',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Utilities',
    ]
)
