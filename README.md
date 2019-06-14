# mmmg_tools

## Install instructions
Optional conda environment creation
```bash
conda create -n mmmg_tools python
source activate mmmg_tools
```
make sure gfortran is default compiler; pymatgen is uncompatible with intel compilers
```bash
pip install -e git+https://github.com/mmm-group/mmmg_tools.git@master#egg=mmmg_tools --src .
```
