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

# Shasta Reservoir: Standalone Simulation

For this first tutorial, we will simulate water temperature in Shasta Reservoir as a single model component - that is, without interacting with any other downstream information or conditions.
Later tutorials will cover a Keswick Reservoir component and methods for coupling these components together.

## A brief introduction to Shasta Reservoir Temperature Simulation

There are many great resources for information on the history and technical details of water temperature management at Shasta Dam and Reservoir -- and how it is commonly simulated in computer models -- so I won't go into much depth here.
For now, the main point to consider is that the temperature of water stored in and released from Shasta Reservoir is affected by several key factors:

  1. The temperature and volume of water flowing into the reservoir
  2. The mixing and stratification (distinctly warmer water separating from cooler water below) of water within the reservoir volume
  3. Heat exchange of water from the area around the reservoir (mainly the surface), and
  4. Selective withdrawal that releases water from different points in the water column to achieve a desired downstream release temperature.

In real life we can see from records of temperature profiles in Shasta Reservoir that each of these factors can have different effects at different times of the year or when the reservoir is empty vs full. 
One key characteristic of Shasta Reservoir is that it statifies seasonally - meaning that a volume of water in the upper part of the vertical profile warms more than water at depth, creating a distinct boundary between two temperature zones. 
The water in this cooler, lower zone is often what is desired for coldwater temperature management for salmon during the hot summer months.

In this version of R3TAM, we use a vertically discretized one-dimensional model we call *ResTemp* (for "Reservoir Temperature") to represent those four factors and their interaction in a reservoir like Shasta.
These factors are addressed by a set of calculations that occur at least once each time step - which for the Shasta *ResTemp* model is one day.
At each time step, the model integrates information read in from user-provided input files, performs calculations, then updates the water temperature and water density conditions (we sometimes refer to these as the system "states").
The output from each time step is a record of those system states as well as other outcomes of interest - primarily the temperature of water released from Shasta Dam (e.g. via penstocks or river outlets) to the tailwater below. 

You can read more about the specifics of the process represenations in the reference sections or the Technical Memorandum document. <<< add link >>>
The construction of the *ResTemp* code - including how the main factors mentioned above are translated to computational terms -- determines the type and form of information and data that need to be provided by a modeler when performing simulations.
We provide a description of these inputs in the sections below on our way toward running a standalone Shasta temperature simulation.

## Model Inputs

The following sections refer to files that support a Shasta Reservoir simulation using the ResTemp model.
You can find these files in the `shasta` folder within the examples provided [on GitHub](https://github.com/james-m-gilbert/r3tam/tree/test-release/examples) or in a zipped folder <<provide-link>>.

### Run Specification File

A ResTemp simulation requires a range of different information and data, including:

  - Reservoir volume, area, elevation relationships
  - Vertical discretization
  - Parameter values
  - Input time series data
  - Release options and characteristics
  - Time period to simulate

This information is provided as input to the model via a in [YAML](https://yaml.org/about/) specification file.
A screenshot of an example YAML specification file for a ResTemp Shasta Reservoir model is shown in {numref}`fig-yamlspec-sha`.
The specification file is structured around *key-value* entries where the *key* is a text entry followed by a colon. 
The *value* entry follows the colon and can take the form of a single value (numeric, text, date, etc) or a path to a separate file that contains other types of information. 

The specification file `shasta_standalone.yaml` provided with the Shasta Reservoir example simulation has all the required information already set, but we provide a brief overview of the types and formats of information below.

```{figure} ./figs/restemp_shasta_specyaml_example.png
:name: fig-yamlspec-sha

Screenshot of top portion of an example ResTemp model specification file in YAML format.
```

#### General Information

The first few entries in the specification file provide general information that applies to the entire simulation and model files.
These are described in {numref}`tbl-restemp-geninput` below:

```{table} ResTemp Model Specification - General Inputs
:name: tbl-restemp-geninput
:align: left

| Entry  | Description        | 
| ---------  | -----------        |
| **RunName**   | Name to be assigned to the model run - output files will use this name, so avoid spaces for better filesystem usability (*required*)|
| **Description** | A longer description of the simulation, helpful to provide notes to yourself about important characteristics of a specific model run (*optional*) |
| **FacilityName** | Name of the reservoir being simulated (*optional*) |
| **SurfaceMethod** | Sets the surface energy exchange method. Currently this should be set to `simple`, but a full energy balance option is being develped (*required - set to `simple`*) |
| **SimulationUnits** | Sets the unit system in which to simulate. Options are `afdegF` for acre-feet and degrees Fahrenheit and `m3degC` for cubic meters and degrees Celsius. For the early release version, the `afdegF` setting is recommended. (*required*) | 
| **Location** | Description, latitute, and longitude provide information about the reservoir. (*optional*)  | 

```

#### Parameters

The second section in the specification file sets the parameter values used by the `ResTemp` simulation components.
The first five (C1-C5) are the primary parameters, and the remaining set of parameters are used to set different options for the simulation.

The parameters and the calibrated values provided as part of the Shasta Reservoir example model are listed in {numref}`tbl-restemp-params`.


```{table} ResTemp Model Specification - Parameters - Values shown are unitless unless otherwise noted
:name: tbl-restemp-params
:align: left

| Parameter | Value | Description |
| --- | --- | --- |
| **C1** | 0.225 | Affects rate of air-water heat exchange at the water surface |
| **C2** | 0.375 | Affects efficiency of heating from absorbed incoming shortwave solar radiation |
| **C3** | 0.980 | Affects efficiency of latent heat exchange due to evaporation |
| **C4** | 0.0525 |Affects rate of heat exchange from inflowing water |
| **C5** | 0.004 | Affects rate of diffusive mixing among layers |
| **DiffuseDist** | 15 (feet) | Distance over which to calculate interlayer diffusion (units set in `DiffuseDistUnits`) |
| **CriticalDepth** | 35 (feet) | Depth from surface at which surface heat exchange goes to zero (units set in `CritcalDepthUnits`) |
| **MeanAlbedo**  | 0.08 | Controls the balance of reflection andabsorption of solar radiation at the water surface |
| **AlbedoAmplitude** | 0.025  | Controls seasonal variability of albedo; if 0, then `MeanAlbedo` applies uniformly for each time step |
| **ExtinctionCoeff** | 0.45 | |
| **SeasSolRad** | 1.0 | An option for modifying solar radiation to account for seasonality; value of 1 means no modification |
| **WindMethod** | 'w2' | Set method by which wind speed is converted to information needed for wind-driven mixing and evaporation; 'w2' uses the same method as in CE-QUAL-W2 |
| **WindAlpha** | 9.45 | Parameter needed for the 'w2' `WindMethod` option |
| **WindBeta** | 0.46 | Parameter needed for the 'w2' `WindMethod` option |
| **WindGamma** | 2.05 | Parameter needed for the 'w2' `WindMethod` option |
| **WindHeight_m** | 10.0 | The height at which input wind speed is measured (assumed in units of meters). This is needed for potential adjustment to 2-m height assumed for wind mixing and evaporation calculations |

```

#### Simulation Properties

The specification file also includes entries to set required configuration information like timing, initial conditions, and simulation mode options.
A screenshot of these entries from the example Shasta specification file is shown in {numref}`fig-timingmode-sha` and explained in more detail below.


 ```{figure} ./figs/restemp_shasta_timingmodeinit_example.png
:name: fig-timingmode-sha

Screenshot of top portion of an example ResTemp model specification file in YAML format.
``` 

##### Timing 

The `*Timing*` entry includes two sub-entries `StartDate` and `EndDate`. 
The input for each should follow a *"month/day/year"* format. 
As an example, to start the simulation on January 1, 2000, one would enter '1/1/2000' as the value to the `StartDate` key.

```{admonition} Some notes on simulation timing
:class: note
:class: dropdown

Check to make sure that the `EndDate` entry is later than the `StartDate`. 
The simulation will not run if this is not the case.

Another thing to keep in mind is that the simulation timing tells the ***ResTemp*** model what time series data to look for as it advances through each timestep. Make sure that input time series provide data for the time range set by the timing.

The ResTemp model assumes a daily time step. The `TimeStepDays` key that is included in example specification files should be set to 1.0 for now. The option for sub-daily time steps is currently under consideration for future development as it requires some substantial changes to the model functions.

```

##### Initial condition 

The *Intial Condition* entry is used to provide input to the model regarding how to intialize the reservoir states (storage volume and temperature). 
Current functionality is limited to reading in and assigning these states from observational data - this is is set by assigning the `Method` sub-entry as *'from_obs'*.
See the <<<< observation data section >>>> for more information about the observation data options and requirements. 

```{admonition} Some notes on simulation timing
:class: note
:class: dropdown

The simulation start date set via `Timing > StartDate` should correspond closely with a date of available observation data - that is, reservoir storage and temperature profiles. In the common case where the desired simulation start date and the available profile data do not exactly align, the model will automatically find the closest date to the set start date and assign the profile data accordingly. This is important to keep in mind when evaluating simulated temperature profiles around these initialization dates.

```

##### Simulation modes


### Reservoir Specifications

The *ResTemp* formulation treats a reservoir as a single volume of water that can vary vertically  but is assumed to be uniform laterally. 
To distinguish this vertical variation in temperature the reservoir volume is assigned a certain number of layers by the modeler. 
The number and thickness (they don't all have to be the same size) of the layers determines how finely or coarsely the resulting temperature simulation will be, and it can affect the accuracy of the simulation.

The text file that specifies these layers for the Shasta example model can be found within the `inputs` folder and is titled `layers42.txt`. 
The file contains a simple listing of index and elevation values (in units of feet) ordered sequentially from low to high. 
Each elevation value sets the bounds of a layer - there are 43 values listed to define the top and bottom elevation of 42 layers.
Note that the layer definitions encompass the elevation range of a full Shasta reservoir - that is, up to 1067 feet above mean sea level. 
This sets the maximum elevation of the reservoir, but the actual volume at any given time step will depend on the mass balance of flows in an out of the reservoir.

The ResTemp model tracks water volume, temperature, and heat exchange at each layer. 
This requires that the volume and area of each layer be defined. 
This is done through two more look-up tables provided in the `inputs` folder:

  - **ShastaElevationAreaTable.txt**: Relates the surface area of each layer top/bottom to the elevation
  - **ShastaElevationCapacityTable.txt**: Sets the total volume of the reservoir when water is filled to that elevation

The ResTemp code reads in these three files (layer specification: `layers42.txt`, elevation-area: `ShastaElevationAreaTable.txt`, & elevation-volume: `ShastaElevationCapacityTable.txt` ) as part of the model initialization and automatically assigns incremental volumes and elevations to each layer that comprises the model. 

### Model Parameters

The ResTemp model uses five main parameters and up to another 10 additional parameters, depending on simulation options, to specify heat exchange and mixing processes. 
These are the parameters that can be adjusted via a calibration process to better fit simulated temperature profiles, fluxes, and release temperatures to observed data.




### Reservoir Inflows





