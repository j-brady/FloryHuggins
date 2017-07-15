# Simple Flory Huggins Phase Diagram Fitting
##

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

`python fit_FT.py`

### Tips



## Theory
```mathjax
 \frac{\Delta F_{mix}}{kT} = \frac{\phi}{N_1} ln \phi + \frac{(1-\phi)}{N_2} ln (1-\phi) + \chi \phi (1-\phi)
```

```mathjax
 \left( \frac{ \partial \Delta F_{mix} }{ \partial \phi } \right) _{\phi=\phi'} = \left( \frac{ \partial \Delta F_{mix} }{ \partial \phi } \right) _{\phi=\phi''}
```
