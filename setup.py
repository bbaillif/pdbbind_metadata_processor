import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pdbbind_metadata_processor",
    version="0.0.1",
    author="Benoit Baillif",
    author_email="benoit.baillif@gmail.com",
    description="Extract metadata from PDBBind index files and give path to files",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bbaillif/pdbbind_metadata_processor",
    project_urls={
        "Bug Tracker": "https://github.com/bbaillif/pdbbind_metadata_processor/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "pdbbind_metadata_processor"},
    packages=setuptools.find_packages(where="pdbbind_metadata_processor"),
    python_requires=">=3.7",
)