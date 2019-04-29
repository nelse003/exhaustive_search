from setuptools import setup
from setuptools import findall

# Not yet working due to libtx.phil import not working?

install_scripts = findall(dir='bin')

setup(
    name="exhaustive_search",
    version="0.1",
    packages=[
        "test",
        "exhaustive",
        "exhaustive.utils",
        "exhaustive.jiffies",
        "exhaustive.plotting",
        "validation",
    ],
    url="",
    scripts= install_scripts,
    license="",
    author="nelse003",
    author_email="nelse003@gmail.com",
    description="Brute force crystallographic occupancy ",
)
