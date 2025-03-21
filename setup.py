from setuptools import setup, find_packages

setup(
    name="ema",
    version='0.1.a1',
    packages=find_packages(),
    install_requires=["pysam>=0.17.0", 
                      "anndata>=0.9.1", 
                      "scanpy>=1.9.3", 
                      "pandas>=1.3.5",
         
                      "scipy>=1.8.0",
                      "igraph>=0.10.6",
                      "louvain>=0.8.0",
                      "pybedtools>=0.9.0"],
                      entry_points={
                          'console_scripts':[
                              'ema=ema.main:main',
                                'ema_switch=ema.switch_test.cli:cli',
                                'ema_merge=ema.merge_bam.cli:cli'
                          ]
                      },
                      author="Amir Amiri Tabat",
                      author_email="amiramiritabat01@gmail.com",
                      url="https://github.com/TRextabat"
                      )
