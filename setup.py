from setuptools import setup, find_packages



setup(

    name='acoutreams',

    version='0.1.0',

    author='Nikita ustimenko',

    author_email='nikita.ustimenko@kit.edu',

    description='T-matrixm scattering code for acoustic computations',

    license = 'MIT',

    long_description=open('README.md').read(),

    long_description_content_type='text/markdown',

    url='https://github.com/NikUstimenko/acoutreams',  

    packages=find_packages(),

    classifiers=[

        'Programming Language :: Python :: 3',

        'License :: OSI Approved :: MIT License',

        'Operating System :: OS Independent',

    ],

    python_requires='>=3.8',

)
