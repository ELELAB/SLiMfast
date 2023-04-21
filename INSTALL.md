# Installing SLiMfast

Last updated: 29/03/2023

## Introduction

`SLiMfast` is a Python package.

You will need Python 3.8 or higher and some open-source Python packages to install `SLiMfast`.

Python is usually available by default in any Linux installation, while the Python packages will be installed automatically during the installation process. Please see the `requirements.txt` file or the `environment.yaml` file for the required Python package list.

## Instructions

Here, we provide instructions for installing `SLiMfast` in a simple Python virtual or `conda` environment.

### Installing in a Python environment with `pip` and `setuptools`

This simple installation procedure requires a Python >= 3.8 distribution to be available. It guides you in installing the `SLiMfast` package in a virtual environment, meaning an instance of Python isolated from your system.

This is not strictly necessary, and `SLiMfast` may be installed system-wide similarly, following steps 4 to 9.

#### Step 1 - Install `virtualenv`

Check if the `virtualenv` Python package is installed in your system.

It is installed if the `virtualenv` command is available. In this case, you can skip this step.

You can usually find it as a package in your distribution if unavailable in your system.

For instance, on Debian-based systems (such as Debian or Ubuntu), this boils down to installing the `python-virtualenv` package.

```shell
sudo apt install python-virtualenv
```

If this is not possible, you may still install the `virtualenv` package for just your local user using `pip`.

```shell
pip install --user virtualenv
```

If the installation is successful, the `virtualenv` command will be available.

#### Step 2 - Create the virtual environment

Create your Python 3.8 virtual environment in a directory of your choice (in this case, it will be `slimfast-env`).

```shell
virtualenv -p /usr/bin/python3.8 slimfast-env
```

You should replace the argument of option `-p` according to the location of your Python installation in the system.

#### Step 3 - Activate the environment

Activate the environment.

```shell
source slimfast-env/bin/activate
```

#### Step 4 - Get SLiMfast

Clone the `SLiMfast` source code from its GitHub repository and enter the local copy of the repository.

```shell
git clone https://github.com/ELELAB/SLiMfast.git
cd SLiMfast
```

If `git` is unavailable, you can download the repository content as a ZIP file from the `SLiMfast` GitHub repository web page.

#### Step 5 - Install the required Python packages

Install all the required Python packages specified in the `requirements.txt` file in the `SLiMfast` directory via `pip`.

```shell
pip install -r requirements.txt
```

#### Step 6 - Install `SLiMfast`

You can now install `SLiMfast`.

```shell
python setup.py install
```

If the installation is successful, the main `SLiMfast` executables (`slimfast_get_structures`, `slimfast_get_homologs`, `slimfast_find_slims`, `slimfast_run_pipeline`) should now be available.

Every time you need to run `SLiMfast` after opening a new shell, just run step 3 beforehand.

### Installing with `conda`

#### Step 1 - Get `SLiMfast`

Clone the `SLiMfast` source code from its GitHub repository.

```shell
git clone https://github.com/ELELAB/SLiMfast.git
```

If `git` is unavailable, you can download the repository content as a ZIP file from the `SLiMfast` GitHub repository web page.

#### Step 2 - Install `conda`

We have provided an `environment.yaml` file that can be used with `conda` to install an environment containing all of `SLiMfast`'s requirements automatically.

To use it, you need to have `conda` installed in your system. Go [here](https://docs.conda.io/en/latest/miniconda.html) for further instructions. Installing the `miniconda` package rather than the full `anaconda` is advised.

Once `conda` is installed on your system, you can create a virtual environment similar to what you would do using the `virtualenv` package, as previously detailed.

#### Step 3 - Create the `conda` environment

Create your `conda` environment starting from the provided environment file.

```shell
conda env create --file SLiMfast/environment.yaml --prefix ./slimfast-env
```

In this case, we ask `conda` to create the environment locally (`--prefix`). This is not strictly necessary.

Note that if your home directory is size-constrained, the installation might fail for lack of space. In that case, you need to select another directory with more space available in which conda will download packages. 

```shell
conda config --add pkgs_dirs /my/other/directory/conda/pkgs
```

You need to change `/my/other/directory/conda/pkgs` to a directory of your choice.

#### Step 4 - Activate the environment

Activate the `conda` environment by running the command line that conda suggests for this purpose at the end of the previous step. It is usually something like this:

```shell
conda activate ./slimfast-env
```

#### Step 5 - Install `SLiMfast`

You can now install `SLiMfast`:

```shell
cd SLiMfast
python setup.py install
```

If the installation is successful, the main `SLiMfast` executables (`slimfast_get_structures`, `slimfast_get_homologs`, `slimfast_find_slims`, `slimfast_run_pipeline`) should now be available.

Every time you need to run `SLiMfast` after opening a new shell, just run step 4 beforehand.