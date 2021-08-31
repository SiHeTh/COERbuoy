from setuptools import setup

setup(
   name='COERbuoy',
   version='0.1.5',
   author='Simon H. Thomas',
   author_email='simon.thomas.2021@mumail.ie',
   packages=['COERbuoy'],
   url='http://coerbuoy.maynoothuniversity.ie',
   license='LICENSE.txt',
   description='A realistic Wave Enegery Converter model to evaluate controllers',
   long_description=open('README.txt').read(),
   install_requires=[
       "numpy",
       "scipy",
       "pandas",
   ],
   include_package_data=True,
)
