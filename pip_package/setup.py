from setuptools import setup

from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name='ThreeDCellComposer',
    version='1.2.1',    
    description='3D cell segmentation by composing 2D segmentations',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/murphylab/3DCellComposer/',
    author='Haoran Chen and Robert F. Murphy',
    author_email='murphy@cmu.edu',
    license='MIT',
    packages=['ThreeDCellComposer'],
    install_requires=['numpy', 'matplotlib', 'scipy', 'skimage', 'deepcell', 'pandas', 'CellSegmentationEvaluator'],

    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',  
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS :: MacOS X',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Image Processing',
    ],
)
