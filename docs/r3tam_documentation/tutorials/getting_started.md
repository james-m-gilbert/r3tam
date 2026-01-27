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

# Prerequisite Knowledge and Skills

The R3TAM framework is written entirely in Python and will require a local (installed on your computer) instance of a Python interpreter to run. 
This tutorial will not cover how to install Python, but there are lots of good resources on the internet that can help you do that. 
I normally use the [Anaconda distribution](https://www.anaconda.com/download) for fairly straightforward Python environment installation.

Although it is not strictly required, it may be helpful to have the `git` version control software installed on your computer. 
If you make changes to the code and want to keep track of them (or want to review the revision history of the existing code base), an understanding of basic version control concepts and how to use `git` commands will be needed.

## A note about Python versions

The R3TAM models were built and tested using Python 3.10.11. 
The code should be compatible with newer versions, but if you want to avoid potential frustrations, I'd suggest sticking with a 3.10.x version.

If you have multiple Python versions on your computer and you'd like to make managing and tracking versions and packages easier, setting up a virtual environment with your base Python version of choice is recommended. 
As an example, I created a virtual environment to test out these tutorials using the following steps (all from within a command prompt/command line window):

  1. Activate Python 3.10.11 using the Conda package manager (`conda activate ~conda-env-name~`)
  2. Navigate to the folder on your computer where you'd like to create the virtual environment, then enter command `python -m venv R3TAMenv` (where R3TAMenv is the name I assigned to the new virtual environment - feel free to change to something different if you prefer)
  3. Once the new environment is created, you should be able to activate it (assuming Windows) using something like Scripts\activate.bat
  4. You should see an indicator in your command window indicating the new environment is active. If not, try de-activating the `conda` environment (`deactivate ~conda-env-name~`) - it may be overriding the new virtual environment.

# Installing the R3TAM Python package

The R3TAM simulation framework is currently distributed via [GitHub](https://github.com/james-m-gilbert/r3tam) and can be installed using the `pip` command that should be part of your Python installation and/or environment.

