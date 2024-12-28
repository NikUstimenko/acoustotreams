from setuptools import setup, find_packages



setup(

    name='acoutreams',

    version='0.1.25',

    author='Nikita Ustimenko',

    author_email='nikita.ustimenko@kit.edu',

    description='A Python package for acoustic scattering based on the T-matrix method',

    license = 'MIT',

    long_description=open('README.md').read(),

    long_description_content_type='text/markdown',

    url='https://github.com/NikUstimenko/acoutreams',

    project_urls={
        "Bug Tracker": "https://github.com/NikUstimenko/acoutreams/issues",
    },  

    packages=find_packages(),

    classifiers=[

        'Programming Language :: Python :: 3',

        'License :: OSI Approved :: MIT License',

        'Operating System :: OS Independent',

    ],

    python_requires='>=3.8',

    install_requires=[
        "numpy",
        "scipy",
        "treams"
    ],

)
