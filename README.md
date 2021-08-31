# COERbuoy
### A realistic Wave Enegery Converter model to evaluate controllers

NOTE: The commands how to run python might differ between systems, and the ones presented here might not work on every machine. Please refer to https://www.python.org/downloads/ to find the correct commands for your system.
<br>
## 1. Repository structure
- COERbuoy (program data)
   - data (WEC data)
   - param (data were generated parameter files are saved)
   - results (folder where results are saved)
   - web (data for GUI)
   - ... (source code)

- examples (examples how to use COERbuoy)
   - custom_WEC (example of a user generated WEC)
   - custom controller (example of a custom controller in python and octave)
   - use_python (how to use COERbuoy from within python)

- dist (files for distribution)

- LICENSE.txt
- README.md (this readme file)
- additional file for the setup configuration

## 2. Installation

Download [python3](https://www.python.org/downloads/).

Install COERbuoy using pip:

| Windows                   |
|:--------------------------|
|`py -m pip install https://github.com/SiHeTh/COERbuoy/raw/main/dist/COERbuoy-0.1.5-py3-none-any.whl`|             

| Linux/MacOS                     |
|:--------------------------------|
| `pip install https://github.com/SiHeTh/COERbuoy/raw/main/dist/COERbuoy-0.1.5-py3-none-any.whl`|

## 3. Run COERbuoy1

Run graphical user interface:

| Windows                   | &nbsp;&nbsp; | Linux/MacOS                     |
|:--------------------------|--------------|:--------------------------------|
|`py -m COERbuoy.GUI`       |              | `python3 -m COERbuoy.GUIServer` |


### 3.1. Command line:

Brettschneider wave with significant wave height 1.5 m and wave energy period of 6 s and a linear generator damping:

| Windows                   |
|:--------------------------|
|`py -m COERbuoy --bretschneider_wave 1.5 6 results.csv linear`|             

| Linux/MacOS                     |
|:--------------------------------|
| `python3 -m COERbuoy --bretschneider_wave 1.5 6 results.csv linear`|

Regular wave with height 1.5 m and period of 6 s and a linear generator damping:

| Windows                   |
|:--------------------------|
|`py -m COERbuoy --regular_wave 1.5 6 results.csv linear`|             

| Linux/MacOS                     |
|:--------------------------------|
| `python3 -m COERbuoy --regular_wave 1.5 6 results.csv linear`|



## 4. Version history

0.1.2 - added compatibility with numpy version > 1.16
