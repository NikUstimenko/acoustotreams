from setuptools import setup, find_packages



setup(

    name='acoustotreams',

    version='0.2.13',

    author='Nikita Ustimenko',

    author_email='nikita.ustimenko@kit.edu',

    description='A Python package for acoustic wave scattering based on the T-matrix method',

    license = 'MIT',

    long_description=open('README.md').read(),

    long_description_content_type='text/markdown',

    url='https://github.com/NikUstimenko/acoustotreams',

    project_urls={
        "Bug Tracker": "https://github.com/NikUstimenko/acoustotreams/issues",
    },  

    packages=find_packages(),

    classifiers=[

        'Programming Language :: Python :: 3.10',

        'Programming Language :: Python :: 3.11',

        'Programming Language :: Python :: 3.12',

        'Programming Language :: Python :: 3.13',

        'License :: OSI Approved :: MIT License',

        'Operating System :: OS Independent',

    ],

    python_requires='>=3.10,<3.14',

    install_requires=[
        "numpy",
        "scipy>=1.14.1,<1.17",
        "treams>=0.4"
    ],

)
