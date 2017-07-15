# Simple Flory Huggins Phase Diagram Fitting
##
Here are the scripts I used to fit the Flory Huggins model to Ddx4 phase diagrams in our submitted manuscript on the structural and hydrodynamic properties of the protein.

Currently, the most simple form of the model is used and salt concentration is not taken into account in the fitting.

## Input data
Tab separated concentration values e.g.

|  Temp  |  dilute  |  std  |  condensed  |  std  |
|--------|----------|-------|-------------|-------|
| deg C  | mg/mL    | mg/mL |  mg/mL      | mg/mL |

where dilute is the concentration of protein in dilute phase, while condensed is the concentration of protein in the condensed (protein rich) phase. Concentrations are in mg/mL.

## Run using a yaml file

create a yaml file called `fit_FH.yml` ...

```yaml
wt:

    model:
        !!python/object:__main__.FloryHuggins
        N1: 153 # e.g. number of amino acids
        N2: 1
        rho: 1400

    data: # you can run over multiple data sets
        - "data_1.txt"
        - "data_2.txt"
        - "data_3.txt"

    p0: [0.0001,0.8] # initial phi1/phi2 params for fitting

    temp_range: [0.,100.] # in degrees C

    outpath: "results_rho1400" # folder where script output is written

```

run `python fit_FT.py`

### Tips



## Theory
![equation](https://latex.codecogs.com/gif.latex?\frac{\Delta&space;F_{mix}}{kT}&space;=&space;\frac{\phi}{N_1}&space;ln&space;\phi&space;&plus;&space;\frac{(1-\phi)}{N_2}&space;ln&space;(1-\phi)&space;&plus;&space;\chi&space;\phi&space;(1-\phi))

![equation](https://latex.codecogs.com/gif.latex?\left(&space;\frac{&space;\partial&space;\Delta&space;F_{mix}&space;}{&space;\partial&space;\phi&space;}&space;\right)&space;_{\phi=\phi'}&space;=&space;\left(&space;\frac{&space;\partial&space;\Delta&space;F_{mix}&space;}{&space;\partial&space;\phi&space;}&space;\right)&space;_{\phi=\phi''})

![equation](https://latex.codecogs.com/gif.latex?\frac{\Delta&space;F_{mix}(\phi')&space;-&space;\Delta&space;F_{mix}(\phi'')}{\phi'&space;-&space;\phi''}&space;=&space;\left(&space;\frac{&space;\partial&space;\Delta&space;F_{mix}&space;}{&space;\partial&space;\phi&space;}&space;\right)&space;_{\phi=\phi'or\phi&space;''})

where

![equation](https://latex.codecogs.com/gif.latex?\left(&space;\frac{&space;\partial&space;\Delta&space;F_{mix}&space;}{&space;\partial&space;\phi&space;}&space;\right)&space;=&space;kT&space;\bigg[&space;\frac{ln\phi}{N_1}&space;&plus;&space;\frac{1}{N_1}&space;-&space;\frac{ln(1-\phi)}{N_2}&space;-&space;\frac{1}{N_2}&space;&plus;&space;\chi(1&space;-&space;2&space;\phi)&space;\bigg])

yielding

![equation](https://latex.codecogs.com/gif.latex?\chi&space;=&space;A&plus;\frac{B}{T}=\frac{\frac{1}{N_1}&space;ln\left(\frac{\phi&space;''}{\phi&space;'}\right)&space;&plus;&space;\frac{1}{N_2}&space;ln\left(\frac{1-\phi&space;'}{1-\phi&space;''}\right)}{2(\phi&space;''&space;-&space;\phi&space;')})
