from setuptools import setup

setup(
   name='COERbuoy',
   version='0.3bN',
   author='Simon H. Thomas',
   author_email='simon.thomas.2021@mumail.ie',
   packages=['COERbuoy'],
   url='http://coer.maynoothuniversity.ie',
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
