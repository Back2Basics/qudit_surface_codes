import setuptools
from pathlib import Path

long_description = Path("README.md").read_text()

setuptools.setup(
    name="Qudit-Surface-Codes-The-Singularity",  # Replace with your own username
    version="0.0.2",
    author="Amelie Schreiber",
    author_email="thesingularity.research@gmail.com",
    description="Surface codes and quantum error correction",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    packages=setuptools.find_packages(),
    install_requires=[
        "matplotlib",
        "networkx",
        "sympy",
        "multiset"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
