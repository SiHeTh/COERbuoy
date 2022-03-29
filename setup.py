from setuptools import setup

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
   name='COERbuoy',
   version='0.3',
   author='Simon H. Thomas',
   author_email='simon.thomas.2021@mumail.ie',
   packages=['COERbuoy'],
   url='http://coerbuoy.maynoothuniversity.ie',
   license='LICENSE.txt',
   description='A realistic Wave Enegery Converter model to evaluate controllers',
   long_description=long_description,
   long_description_content_type="text/markdown",
   install_requires=[
       "numpy",
       "scipy",
       "pandas",
   ],
   include_package_data=True,
)
