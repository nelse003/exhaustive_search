from setuptools import setup
from setuptools import findall
from setuptools import find_packages

# Not yet working due to libtx.phil import not working?

setup(
    name="exhaustive",
    version="0.1",
    packages=find_packages(),
    url="",
    scripts=findall(dir="bin"),
    license="",
    author="nelse003",
    author_email="nelse003@gmail.com",
    description="Brute force crystallographic occupancy ",
)
