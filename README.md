# COERbuoy
A realistic Wave Enegery Converter model to evaluate controllers

# Installation
Download python3: https://www.python.org/downloads/

Install COERbuoy using pip: pip install https://github.com/SiHeTh/COERbuoy/raw/main/dist/COERbuoy-0.1.0-py3-none-any.whl

# Run COERbuoy1
## Bretschneider sea state:

With significant wave height 1.5 and wave energy period of 6 and a linear generator damping:

Windows: python -m COERbuoy --bretschneider_wave 1.5 6 results.csv linear

Linux/MacOS: python3 -m COERbuoy --bretschneider_wave 1.5 6 results.csv linear


## Regular wave:

With wave height 1.5 and wave period of 6 and a linear generator damping:

Windows: python -m COERbuoy --regular_wave 1.5 6 results.csv linear

Linux/MacOS: python3 -m COERbuoy --regular_wave 1.5 6 results.csv linear
