import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="olivar",
    version="1.2.1",
    author="Michael X. Wang",
    author_email="xw66@rice.edu",
    description="Olivar PCR tiling design",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://gitlab.com/treangenlab/olivar",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL-3.0 license",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=['olivar'],
    python_requires=">=3.8",
)
