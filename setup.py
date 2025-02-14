from setuptools import setup, find_packages

setup(
    name="needlemann_wunsch", 
    version="0.1.0",
    author="Justin Sim",
    author_email="justin.sim@ucsf.edu",
    description="Implementation of the Needleman-Wunsch algorithm for global sequence alignment",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/justinsim12/HW5-NW",  
    install_requires=[
        "numpy"
    ],
    extras_require={
        "dev": ["pytest"],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)