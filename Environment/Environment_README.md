# Environment Creation

Prerequisite: miniconda3 installed in your system.
If not, use the link: https://www.anaconda.com/download/success?reg=skipped-miniconda
and install it in your system. 

Use the installed specific to your system OS: Windows or Linux

## For Windows OS:

### Installation: 
1) Use the file "pcb_v1.yml" to create the conda environment in your system.
2) Use the following command to install the environment:

```bash
conda env create -f pcb_v1.yml
```
#### NOTE:  You need to be in the same directory while installing the env, 
otherwise give the entire path to the file like:

```bash
conda env create -f c:/to/the/path/pcb_v1.yml
```
Your environment will be installed.

### Activating conda env:

Use the following command to activate the pcb_v1 environment:

```bash
conda activate pcb_v1
```
Use the following command to deactivate the pcb_v1 environment:

```bash
conda deactivate
```
