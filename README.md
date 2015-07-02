# MicrobiomeSim

This repository contains a set of programs and scripts for simulating microbial communities under 
known perturbations to abundances of individuals taxa and changes in overall diversity. It allows 
us to determine if a change of given magnitude is likely to be detected. 

##GenerateSimVar.py

This program takes as input an file in csv format of taxa and their frequencies across samples. This 
may be obtained from a microbiomics 16S rRNA study for example. Using this to define the control 
community it then generates simulated data sets of control and perturbed data. 

The constants for generating the simulated data are hard coded in the file.

```
class Constants(object):
    LOG_SCALE_POP = 0.0
    SAMPLE_SIZE   = 2000
    N_CONTROL     = 50
    N_TREATMENT   = 50
    N_CHANGE      = 25
    THRESHOLD     = -8.75
    MAX_CHANGE    = np.log(10.0)
    MIN_CHANGE    = -np.log(10.0)
```




```
python ./src/GenerateSimVar.py data/hmpv13_R_crease_PPS_R1e3.csv data/hmpv13_R_crease_PPS_R1e3_test
```
