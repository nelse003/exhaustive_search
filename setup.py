from setuptools import setup

# Not yet working due to libtx.phil import not working?

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
    entry_points={
        "console_scripts": [
            "exhaustive=exhaustive.exhaustive.__main__:main",
            "exhaustive_multiple_sample=exhaustive.exhaustive_multiple_sampling.__main__:main",
        ]
    },
    license="",
    author="nelse003",
    author_email="nelse003@gmail.com",
    description="Brute force crystallographic occupancy ",
)
