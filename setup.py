import setuptools

long_description = "a scanpy-based single-cell crosstalk analysis package"

setuptools.setup(
  name="depair",
  version="0.0.1",
  author="Sijie Chen",
  author_email="chansigit@gmail.com",
  description="A scanpy-based single-cell crosstalk analysis package",
  long_description=long_description,
  long_description_content_type="text/markdown",
  url="https://github.com/chansigit/depair",
  packages=setuptools.find_packages(),
  classifiers=[
  "Programming Language :: Python :: 3",
  "License :: OSI Approved :: MIT License",
  "Operating System :: OS Independent",
  ],
)