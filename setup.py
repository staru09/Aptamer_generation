from setuptools import setup, find_packages

setup(
    name="aptamer_generator",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "biopython",
        "matplotlib",
    ],
    author="Aru Sharma",
    author_email="arusharmazxx000@gmail.com",
    description="A package for generating candidate aptamers",
    keywords="aptamer, bioinformatics, DNA",
    url="https://github.com/staru09/aptamer-generator",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
    ],
    python_requires=">=3.8",
    package_data={
        "aptamer_generator": ["examples/*.py"],
    },
)