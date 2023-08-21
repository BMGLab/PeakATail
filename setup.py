from setuptools import setup, find_packages

setup(
    name="ema",
    version='0.1.a1',
    packages=find_packages(),
    install_requires=["pysam>=0.17.0", 
                      "anndata>=0.9.1", 
                      "scanpy>=1.9.3", 
                      "pandas>=1.3.5", 
                      "scipy>=1.8.0"],
                      author="Amir Amiri Tabat",
                      author_email="amiramiritabat01@gmail.com",
                      url="https://github.com/TRextabat"
                      )
