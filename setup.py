from setuptools import setup

setup(
    name='bn256',
    version='0.1.0',
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    packages=['bn256'],
    url='https://github.com/hanzhichao/bn256.git',
    license='MIT License',
    author='Han Zhichao',
    author_email='superhin@126.com',
    description='python bn256',
    install_requires=[],
    keywords=["bn256", "python-bn256", "hibe"],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Topic :: Software Development :: User Interfaces',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ]
)
