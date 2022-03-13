# COERbuoy
### A realistic Wave Enegery Converter model to evaluate controllers

#### Fore more information, have a look at the [manual](https://github.com/SiHeTh/COERbuoy/raw/main/manual.pdf) ####
#### Or see the basic_usage jupyter notebook. ####

##### Version 0.3 beta #####

NOTE: The commands how to run python might differ between systems, and the ones presented here might not work on every machine. Please refer to https://www.python.org/ to find the correct commands for your system.
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

| Variant 1          |
|:--------------------------|
|`py -m pip install COERbuoy`|             

| Variant 1              |
|:--------------------------------|
| `python3 -m pip install COERbuoy`|

## 3. Run COERbuoy

Run graphical user interface:

| Variant 1                    | &nbsp;&nbsp; | Variant 1              |
|:--------------------------|--------------|:--------------------------------|
|`py -m COERbuoy.GUI`       |              | `python3 -m COERbuoy.GUI      ` |


### 3.1. Command line:

Brettschneider wave with significant wave height 1.5 m and wave energy period of 6 s and a linear generator damping:

| Variant 1                    |
|:--------------------------|
|`py -m COERbuoy --bretschneider_wave 1.5 6 results.csv linear`|             

| Variant 1             |
|:--------------------------------|
| `python3 -m COERbuoy --bretschneider_wave 1.5 6 results.csv linear`|

Regular wave with height 1.5 m and period of 6 s and a linear generator damping:

| Variant 1                    |
|:--------------------------|
|`py -m COERbuoy --regular_wave 1.5 6 results.csv linear`|             

| Variant 1              |
|:--------------------------------|
| `python3 -m COERbuoy --regular_wave 1.5 6 results.csv linear`|



## 4. Version history
- 0\.1\.2 
   - added compatibility with numpy version > 1.16
- 0\.2
   - manual added
   - settings tab added in GUI
   - simplify use of custom controller (via new folder and automatical detection of ocatve/python interpreter)
