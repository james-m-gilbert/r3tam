---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---
# Getting Started with the Rapid Reservoir-River Temperature Assessment Modeling (R3TAM) Framework



## Prerequisite Knowledge and Skills

The R3TAM framework is written entirely in Python and will require a local (installed on your computer) instance of a Python interpreter to run. 
This tutorial will not cover how to install Python, but there are lots of good resources on the internet that can help you do that. 
I normally use the [Anaconda distribution](https://www.anaconda.com/download) for fairly straightforward Python environment installation.

Although it is not strictly required, it may be helpful to have the `git` version control software installed on your computer. 
If you make changes to the code and want to keep track of them (or want to review the revision history of the existing code base), an understanding of basic version control concepts and how to use `git` commands will be needed.

```{note}
##### Python Versions
The R3TAM models were built and tested using Python 3.10.11. 
The code should be compatible with newer versions, but if you want to avoid potential frustrations, I'd suggest sticking with a 3.10.x version.

If you have multiple Python versions on your computer and you'd like to make managing and tracking versions and packages easier, setting up a virtual environment with your base Python version of choice is recommended. 


  ```

There is no graphical user interface (GUI) to accompany the models, so all interactions with the model are through code, configuration, and data files. 
The tutorials, for the most part, assume that the use has some experience and comfort with common tasks like navigating via a command prompt and editing/running code in an interpreter. 

## Setting up a Python Environement

If you work with Python frequently and have a mixed set of versions and environments, it's cleaner (and can save you some future frustrations) to create a "virtual environment" specifically for the task at hand. 
In this case, that would mean creating an virtual environment for the R3TAM modeling of Shasta, Keswick, and the Sacramento River.

There are good online resources that can give you an introduction to Python virtual environments (e.g. [Creation of Virtual Environments - Python 3.10](https://docs.python.org/3.10/library/venv.html)). Check these out if you want to learn more.

As an example, here are steps I used to create a virtual environment for running R3TAM simulations (all from within a command prompt/command line/console window):

  1. Activate Python 3.10.11 using the Conda package manager (`conda activate ~conda-env-name~`, where *~conda-env-name~* is the name of your base Conda environment, `py310`, for example)
  2. Navigate to the folder on your computer where you'd like to create the virtual environment, then enter command `python -m venv R3TAMenv` (where *R3TAMenv* is the name I assigned to the new virtual environment - feel free to change to something different if you prefer)
  3. Once the new environment is created, you should be able to activate it (assuming Windows) using something like Scripts\activate.bat
  4. You should see an indicator in your command window indicating the new environment is active. If not, try de-activating the `conda` environment (`deactivate ~conda-env-name~`) - it may be overriding the new virtual environment.


An example of what these steps look like in a command prompt window is shown below:

```{code} console

conda activate py310

cd Modeling\WaterTemperature\

python -m venv R3TAMenv

```

Note that the second command will create the virtual environment in a folder called `R3TAMenv` in the current working directory of the command console. 
If you want to put the environment in a specific location (say, in a folder where you will keep all your temperature model runs), you can either navigate to that directory or include the full path specifying the intended location of `R3TAMenv`.

Once you have the virtual environment created, you can activate it so that any Python commands you enter (or scripts you run), will reference this new Python environemnt. 
The activation commands vary based on your platform:
 
  - Windows: (cmd.exe) `C: <venv>\Scripts\activate.bat`
  - Linux/macOS: (bash): `source <venv>/Scripts/activate`

One final step - I recommend setting up a IPython/Jupyter kernel so you can interact with the R3TAM code in an IDE environment. To do that, you'll need to install additional libraries. I tend to use the Spyder IDE most frequently, so the command to install necessary components is: `pip install spyder-kernels==3.0.*`
  
Once you've successfully activated your R3TAM virtual environment and isntalled the `spyder-kernels`, you are ready for installation!


# Installing the R3TAM Python package

The R3TAM simulation framework is currently distributed via [GitHub](https://github.com/james-m-gilbert/r3tam) and can be installed using the `pip` command that should be part of your Python installation and/or environment.

You can install the code into your selected local environment using the following command:

```{code} console
:name: github-install-r3tam

pip install git+https://github.com/james-m-gilbert/r3tam@test-release
```

The first part of the command (`pip`) calls the Python package manager and, as long as your R3TAM environment is active, the `pip` commands will apply just to that environment.
The `install` command indicates a package is to be installed.
The third part of the command specifies that we're installing from a GitHub repository. In this case we're installing a specific branch or version - the "test-release" branch, as indicated by that name following the *@* after `r3tam`.

On running this command in the command line, you should see a bunch of test scroll by - downloading, installing, and hopefully a "Successfully built r3tam" followed by "Successfully installed ..." a listing of the dependencies on which the R3TAM code relies. 
If you don't see any error messages, the package should be installed and ready to test out. 
See the tutorials for Shasta Reservoir, Keswick Reservoir, Sacramento River, and coupled simulation to learn how to run the different model components.