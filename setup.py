from setuptools import setup, find_packages

setup(
    packages=find_packages(include=["gget", "gget.*"]),
    include_package_data=True,
    entry_points={
        "console_scripts": ["gget=gget.main:main"],
    },
)
