from setuptools import setup

setup(
   name='COERbuoy',
   version='0.1.0',
   author='Simon H. Thomas',
   author_email='simon.thomas.2021@mumail.ie',
   packages=['COERbuoy'],
   #scripts=['bin/script1','bin/script2'],
   url='coerbuoy.maynoothuniversity.ie',
   license='LICENSE.txt',
   description='A realistic Wave Enegery Converter model to evaluate controllers',
   long_description=open('README.txt').read(),
   install_requires=[
       "numpy",
       "scipy",
       "pandas",
   ],
   include_package_data=True,
   package_data={'':['/data/*']},
)