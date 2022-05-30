# Snakemake for reproducible analyses

This repository is used as a base for the workshop "Introduction to Snakemake for reproducible analyses" organized at SIB-Days on June 15th, 2022. It is an updated version of a previous workshop organized at [workshop] organized at [PyPharma 2019](https://github.com/RomainFeron/workshop-snakemake-pypharma2019), [SIB-days 2020](https://github.com/RomainFeron/workshop-snakemake-sibdays2020), and at [UNIL in October 2021](https://github.com/RomainFeron/workshop-snakemake-unil2021). The repository contains:

- A conda environment with all the dependencies required during the workshop
- A reference implementation of solutions for all exercises
- A 'workflow' directory with sample data in which participants will implement the exercises

The structure of the workshop loosely follows that of the official [Snakemake tutorial](https://snakemake.readthedocs.io/en/v7.8.0/tutorial/tutorial.html), with a few modifications.

All the exercises and material required to complete them are available on the repository's [wiki](https://github.com/RomainFeron/workshop-snakemake-sibdays2022/wiki).

## Setup the workshop environment

## Note for Windows user

Although Snakemake and Conda can work on Windows, the shell and general environment are very different from unix-based operating systems. If you are using a Windows machine for this workshop, please refer to the [Windows installation instructions](https://snakemake.readthedocs.io/en/v6.9.1/tutorial/setup.html#setup-on-windows) in the official documentation. We recommend you to setup the WSL if you can, as it is the most efficient way to run Snakemake on Windows, and it will be useful for many other Bioinfomatics applications in the future.

### Installing Conda

Detailed instructions [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

**Download the installer script and run it:**

*Linux:*

```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

*MacOS:*

```bash
curl https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh
```

**Follow instructions from the prompt:**

- Accept the License agreement (`yes`)
- Chose the installation path (default is `/home/<user>/miniconda3`)
- Wait for the packages to be installed (this can take some time)
- Chose whether the installer should initialize conda for your shell (recommended: `yes`)

**Restart your shell:**

```bash
source ~/.bashrc
```

**Note:** If you're not using bash as your shell, source the appropriate file (*e.g.* `~/.zshrc` if you're using zsh)

To verify that the installation process completed correctly, run:

```bash
conda env list
```

The output should look like this (username and hostname will be different):

```bash
(base) user@host:~$ conda env list
# conda environments:
#
base                  *  /home/nbuser/miniconda3
```

**Update conda:**

The Conda version from the official installer is not always the latest update. To make sure Conda is up-to-date, run:

```bash
conda update conda
```

**Install mamba:**

Mamba is a drop-in replacement of the conda executable which is much more efficient to resolve dependencies and install packages. It is used by Snakemake by default and we recommend using it to install packages in the future. Install it with:

```bash
conda install -c conda-forge mamba
```

After the install process is done, you can simply replace `conda` with `mamba` in most conda commands. There's no real point to do that for fast commands like creating an empty environment or listing environments, but it is much faster to use mamba when installing packages; in the commands below, we will use mamba when appropriate.

### Clone the workshop's repository

**With SSH (recommended in general for GitHub, need to set it up if you haven't already):**

```bash
git clone git@github.com:RomainFeron/workshop-snakemake-sibdays2022.git
```

**With HTTPS (default, no setup required):**

```bash
git clone https://github.com/RomainFeron/workshop-snakemake-sibdays2022.git
```

### Create a Conda environment with all software required by the workshop

We provide an environment file `workshop.yaml` that contains all software required to complete the workshop.

Navigate to the workshop's base directory. If you followed the previous instructions exactly, you can do with:

```bash
cd workshop-snakemake-sibdays2022
```

Create the conda environment from the environment file with mamba:

```bash
mamba env create -f workshop.yaml
```

This step usually takes some time (up to 20-30 minutes), as there is a lot of dependencies to install.

Once the environment is created, activate it with:

```bash
conda activate snakemake-workshop
```

You can now run Snakemake and complete the workshop's exercises.

### Additional note: creating an Conda environment for Snakemake

For this workshop, for the sake of simplicity, we created a single Conda environment that contains both Snakemake and all the tools used in the workflow. In practice, this is not the correct way to use Conda with Snakemake workflows, as we will discuss during the workshop. This step is not part of the setup process, but if you want to use Snakemake by yourself in the future, the recommended way to run it is to create a conda environment specifically and only for Snakemake:

```bash
# Create a new empty environment called "snakemake"
conda create --name snakemake
# Activate the environment "snakemake"
conda activate snakemake
# Install snakemake from the Bioconda channel (conda-forge contains dependencies)
mamba install -c conda-forge -c bioconda snakemake
```

You can now activate the environment snakemake and run Snakemake from it. It is advised to keep the environment as clean as possible, *i.e.* only install software related to running snakemake in general, not software specifically for your workflow.
