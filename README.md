# Rodriguez Lab Notebook

This repository stores the code Franklin Hiciano used to generate data & results in the [Rodriguez Lab](https://oscarlr.github.io/), which aims to understand the IG locus. The code is auto-generated from a Colab notebook using custom scripts.

## Installation
```
pip install --no-cache-dir --force-reinstall "git+https://github.com/franklin-hiciano/rodriguez-lab-notebook.git@main"
```
In Jupyter/Colab notebooks it is recommended to use:
```
!{sys.executable} -m pip install --no-cache-dir --force-reinstall "git+https://github.com/franklin-hiciano/rodriguez-lab-notebook.git@main"
```

## Usage
Ex:
```python
>>> from FranklinLabNB.data.SabDabRetriever import SabDabRetriever
```
```python
>>> s = SabDabRetriever()
>>> summary_file = s.get_summary_file()
Saved SabDab summary file to sabdab_summary_2025-08-22_04-07-36.tsv.
>>> structures_dir = s.get_structures()
Saved SabDab structures to sabdab_complexes_2025-08-22_03-53-15
```
If the instance is deleted you can use manual flags, like:
```python
structures_dir = s.get_structures(summary_file="sabdab_summary_2025-08-22_04-07-36.tsv")
```
or
```python
structures_dir = s.get_structures(summary_file="sabdab_summary_2025-08-22_04-07-36.tsv")
```

## Authors
Franklin Hiciano
[fhiciano5@gmail.com]

## Acknowledgments

1. CEYE, who made this all possible
2. Beloved lab members
   - Oscar L. Rodriguez, PI
   - Jackie Woodruff
   - Igor Lobo
   - Jacquelyn Willis
   - Nico Araya
   - Nefte Mendoza
3. Icahn School of Medicine at Mount Sinai
