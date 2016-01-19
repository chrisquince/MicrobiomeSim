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
    SAMPLE_SIZE   = 2000
    N_CONTROL     = 50
    N_TREATMENT   = 50
    N_CHANGE      = 25
    THRESHOLD     = -8.75
    MAX_CHANGE    = np.log(10.0)
    MIN_CHANGE    = -np.log(10.0)
```

These define the following:

* SAMPLE_SIZE - the size of each simulated data set
* N_CONTROL - the number of control samples
* N_TREATMENT - the number of treatment samples
* N_CHANGE - number of OTUs chosen randomly to be perturbed
* THRESHOLD - minimum log abundance for perturbed OTU
* MAX_CHANGE - maximum change in OTU mean abundance
* MIN_CHANGE - minimum change in OTU mean abundance

The simulated data sets are generated by first fitting a log-normal to each OTU abundance across sites. These are used for the control samples. For the 
perturbed samples a random factor is added to each OTU mean abundance. The samples are then generated by multinomial sampling.

The script is run by simply specifying the empirical input file and the output file stub to be generated.

```
python ./src/GenerateSimVar.py data/hmpv13_R_crease_PPS_R1e3.csv data/hmpv13_R_crease_PPS_R1e3_test
```

##GenerateDPerturb.py

This program takes as input an file in csv format of taxa and their frequencies across samples. This 
may be obtained from a microbiomics 16S rRNA study for example. Using this to define the control 
community it then generates simulated data sets of control and perturbed data. In this case perturbation 
corresponds to a reduction in diversity prior to subsampling to randomly selected sizes. The constants for 
generating the simulated data are command line arguments of the program.

```
usage: GenerateDPerturb.py [-h] [-s RANDOM_SEED] [-d DIVERSITY_REDUCTION]
                           [-ca CONTROL_A] [-ta TREATMENT_A] [-cb CONTROL_B]
                           [-tb TREATMENT_B] [-l BASE_LEVEL]
                           input_file output_file_stub
```

These define the following:
parser.add_argument('-d','--diversity_reduction',default=0.5, type=float,
        help=("fraction of OTUs removed by perturbation"))

    parser.add_argument('-ca','--control_a',default=4.0, type=float,
        help=("shape paramter for control gamma distn"))

    parser.add_argument('-ta','--treatment_a',default=4.0, type=float,
        help=("shape paramter for treatment gamma distn"))

    parser.add_argument('-cb','--control_b',default=2000, type=float,
        help=("scaleXshape paramter for control gamma distn"))

    parser.add_argument('-tb','--treatment_b',default=2000, type=float,
        help=("scaleXshape paramter for control gamma distn"))        

    parser.add_argument('-l','--base_level',default=10000, type=int,
        help=("specifies initial sub sampling"))

* DIVERSITY_REDUCTION - the fraction of OTUs to randomly selected and removed - default 0.5
* CONTROL_B - beta/alpha i.e. mean parameter for control sample sizes - default 2000
* TREATMENT_B - beta/alpha i.e mean parameter for treatment sample sizes - default 2000
* CONTROL_A - alpha parameter i.e. std. dev. = mean/sqrt(alpha) for control sample sizes - default 4
* TREATMENT_A - alpha parameter i.e. std. dev. = mean/sqrt(alpha) for treatment sample sizes - default 4

The simulated data sets are generated by first generating sample sizes for the control and treatment 
samples from a Gamma distribution. Each control sample is then generated by subsampling each empirical 
sample to that size. The same is done for the treatment samples except in that case a fraction of the 
OTUs are first removed. The samples are then generated by multinomial sampling.

The script is run by simply specifying the empirical input file and the output file stub to be generated.

```
python ./src/GenerateDPerturb.py data/s10e3_R0.csv data/s10e3_R0_out
```
