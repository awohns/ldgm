from setuptools import setup

with open("README.md") as fh:
    long_description = fh.read()

setup(
    name="ld_graph",
    author="Anthony Wilder Wohns",
    author_email="awohns@gmail.com",
    description="Create a graphical model of SNP conditional dependence",
    long_description=long_description,
    packages=["ld_graph"],
    long_description_content_type="text/markdown",
    python_requires=">=3.4",
    entry_points={
        "console_scripts": [
            "ld_graph=ld_graph.__main__:main",
        ]
    },
    setup_requires=["setuptools_scm"],
    install_requires=[
        "tskit>=0.3.0",
        "flake8",
        "numpy",
        "networkx",
        "pytest",
        "pandas"
    ],
    project_urls={
        "Source": "https://github.com/awohns/ld_graph",
        "Bug Reports": "https://github.com/awohns/ld_graph/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3.5",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
    ],
    use_scm_version={"write_to": "ld_graph/_version.py"}
)
