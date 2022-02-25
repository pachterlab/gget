from setuptools import setup, find_packages

def read(path):
    with open(path, 'r') as f:
        return f.read()

long_description = read('README.md')

setup(
    name='gget',
    version='0.0.5',
    license='MIT',
    author='Laura Luebbert',
    author_email='lauraluebbert@caltech.edu',
    maintainer='Laura Luebbert',
    maintainer_email='lauraluebbert@caltech.edu',
    description='Query Ensembl for genes using free form search words, look up genes/transcripts by Ensembl ID or fetch the latest FTPs by species.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=find_packages(),
    install_requires=[
        'mysql-connector-python',
        'beautifulsoup4>=4.10.0',
      ],
    url='https://github.com/lauraluebbert/gget',
    keywords='gget',
    entry_points={
        'console_scripts': ['gget=gget.gget:main'],
    },
    classifiers=[
        'Development Status :: 5 - Production/Stable',
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
