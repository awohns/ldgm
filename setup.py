from setuptools import setup

with open("README.md") as fh:
    long_description = fh.read()

setup(
    name="ldgm",
    author="Anthony Wilder Wohns, Luke O'Connor, Pouria Salehi Nowbandegani",
    author_email="awohns@gmail.com",
    description=" Sparse and accurate modeling of linkage disequilibrium",
    long_description=long_description,
    packages=["ldgm"],
    long_description_content_type="text/markdown",
    python_requires=">=3.4",
    entry_points={
        "console_scripts": [
            "ldgm=ldgm.__main__:main",
        ]
    },
    setup_requires=["setuptools_scm"],
    install_requires=[
        "tskit>=0.3.0",
        "flake8",
        "numpy",
        "networkx",
        "pytest",
        "pandas",
    ],
    project_urls={
        "Source": "https://github.com/awohns/ldgm",
        "Bug Reports": "https://github.com/awohns/ldgm/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3.5",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
    ],
    use_scm_version={"write_to": "ldgm/_version.py"},
)
