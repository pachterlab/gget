from setuptools import setup, find_packages

def read(path):
    with open(path, 'r') as f:
        return f.read()

long_description = read('README.md')

setup(
    name='gget',
    version='0.0.4',
    license='MIT',
    author='Laura Luebbert',
    author_email='lauraluebbert@caltech.edu',
    description='Query Ensembl for genes using free form search words (gget search) or fetch FTP download links by species (gget FetchTP).',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=find_packages(),
    install_requires=[
          'mysql-connector-python',
      ],
    url='https://github.com/lauraluebbert/gget',
    keywords='gget',
    entry_points={
        'console_scripts': ['gget=gget.gget:main'],
    },
)
