from setuptools import setup, find_packages


with open("README.md", "r") as f:
    long_description = f.read()


with open("gget/__init__.py", "r") as f:
    init = f.readlines()

for line in init:
    if '__author__' in line:
        __author__ = line.split("'")[-2]
    if '__email__' in line:
        __email__ = line.split("'")[-2]
    if '__version__' in line:
        __version__ = line.split("'")[-2]

setup(
    name='gget',
    version=__version__,
    author=__author__,
    author_email=__email__,
    description='Query Ensembl for genes using free form search words (gget search) or fetch FTP download links by species (gget FetchTP).',
    long_description=long_description,
    long_description_content_type='ext/markdown',
    packages=find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ),
)
