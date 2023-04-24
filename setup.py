import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="analysis_functions",
    author="Christina Asmus",
    author_email="christina.asmus@hereon.de",
    description="irrigation functions",
    python_requires=">=3.7",
    package_dir={"analysis_functions": "analysis_functions"}
    # install_requires=['Pillow'],
)
