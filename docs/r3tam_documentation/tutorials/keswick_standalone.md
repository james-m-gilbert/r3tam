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

# Tutorial: Simulating Keswick Reservoir - an example fully-mixed tailbay reservoir type

The three main component models in the R3TAM framework represent different conceptualizations of how water temperature changes through managed resevoir-river systems.
Tailbay or re-regulating reservoirs tend to be smaller volume impoundments that sit below a larger reservoir and serve the purpose of smoothing out variations in releases.
They may also serve to collect and manage inflows from multiple sources, either natural or managed. 
For purposes of modeling, such reservoirs are assumed to mix fully along the major axes of flow, resulting in released water temperatures that match the fully-mixed condition in the reservoir.
It is also assumed in this conceptualization that there is no significant vertical or lateral variation in temperature and that there is no selective withdrawal capability at the dam (after all, without vertical stratification, there is no water of different temperature to selectively withdraw).

Keswick Reservoir, downstream of Shasta Reservoir on the Sacramento River in northern California, is a good example of this tailbay-type reservoir. This tutorial walks through the setup and execution of an R3TAM ~LongTemp~ mixed model of Keswick Reservoir for the recent historical period.


## Keswick Reservoir Specifications

asdfas. 


## Model Inputs

asdf


### Meteorology

asdf

### Parameters

asdf



## Model Setup

Up to this point we have not run any code - we've just been setting up the input files that a simulation will use. 
In order to run our example Keswick Reservoir model, we'll need to write a Python file (*.py) that connects the R3TAM framework models to the input data we've created. 
The following steps assume you have already installed the R3TAM code in a Python environment on your computer, and that you've activated that environment prior to running any of this code (or are referencing the correct environment kernel) in your IDE or editor of choice.

```{note}
The following guide and information assumes the code is being executed in a script that is executed as a single file - either from a command line or from within a code editor (e.g. Spyder)
The code could just as easily be transferred to a Jupyter notebook environment and executed by cell to achieve the same results.
```

### Import Module 

We will be using the `longtemp` component of R3TAM for the Keswick reservoir model, so the first step in setting up a run is to import this component from the package.
If you were doing this from scratch, you'd first, open (or create) a new Python (\*.py) file in your favorite editor or IDE  and name the file something descriptive like *keswick_tutorial.py*.
We've already created this file - it's available in the **examples/keswick** directory of the R3TAM [GitHub repository](https://github.com/james-m-gilbert/r3tam) if you want to follow along - otherwise you can build the script according to the steps described below.
Assuming you've successfully installed the R3TAM package in your chosen Python environment, you can import the `longtemp` module using the following line:

```python
from r3tam.longtemp import LongTemp as lt
```

This import syntax allows us to use `lt` as shorthand for the full `longtemp` module and `LongTemp` class - now we can reference functions related to this model component with a bit less typing.
The first function associated with the `lt` object that we'll use here is .`initialize_longmod`.
This function takes as an input the path to the *.yaml configuration file we created and, after successfully reading in all of that information, returns a `LongTemp` model object that we can then use for making simulations.
Assuming we're using the configuration file developed above (and provided in the **examples/keswick** directory), and that our interpreter is set to the **examples\keswick** folder as the current working directory, we can initialize the Keswick model using the following commands:

```python
input_fp = r"kwk_input.v20240729.yaml"

kwk = lt.initialize_longmod(input_fp, show_init=True)
```

```{margin}
  {hint} 
  If you have issues with the `initialize_longmod` function *not* finding the configuration file, try providing the full file path.
```

The first line sets the path to the input specification file.
The next line calls the initialization function and creates our Keswick model object -- in this case we assign that to the variable `kwk`.


