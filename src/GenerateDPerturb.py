import sys, getopt
import pandas as p
import numpy as np
import random as rnd
import argparse

class Constants(object):
    DIVERSITY_REDUCTION   = 0.5
    #CONTROL_B             = 5000
    CONTROL_B             = 2500
    TREATMENT_B           = 2500
    #TREATMENT_B           = 5000
    CONTROL_A              = 6.25
    #CONTROL_A             = 25
#    TREATMENT_A           = 25 
    TREATMENT_A           = 6.25

#generate multinomial samples from matrix of relative abundances
def generate_MN(rel_abund_matrix,n_Multi_Samples):

    ra_dim = rel_abund_matrix.shape
    nSamples = ra_dim[0]
    nOTUs = ra_dim[1]

    otus = np.zeros(shape=(nSamples,nOTUs),dtype=np.int)
    for i in range(nSamples):
        otus[i,:] = np.random.multinomial(n_Multi_Samples[i],rel_abund_matrix[i,:])

    return otus

def set_random_state(seed):
    ERROR="'{0}' should be converatable to integer".format(seed)
    try:
        seed = int(seed)
        if seed < 0:
            raise ArgumentTypeError("'" + seed + "' should be >= 0")
        elif seed == 0:
            seed = randint(2,10000)
        return seed
    except ValueError as e:
        raise ArgumentTypeError(ERROR)

def generate_sample_df(control_otus,otu_labels,sample_labels):
    sampled_otus = np.vstack((control_otus))

    sampled_otus_df = p.DataFrame(sampled_otus,columns=otu_labels)

    sampled_otus_df['Samples'] = sample_labels

    cols = sampled_otus_df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    sampled_otus_df = sampled_otus_df[cols]
    return sampled_otus_df

def generate_logSS(a,b,nSamples):
    
    Z = np.random.gamma(a,b,size=(nSamples))

    ZV = Z.astype(int)

    return ZV

def perturb_diversity(input_otus,input_otus_P,nSamples,nOTUs):

    perturb_otus_P = np.copy(input_otus_P)
    perturb_otus   = np.copy(input_otus)
    for s in range(nSamples):
        #get list of nonzero indices
        nonzero_s = np.nonzero(input_otus_P[s,:])
        nP = len(nonzero_s[0])
        nS = int(Constants.DIVERSITY_REDUCTION*nP)
        remove = rnd.sample(nonzero_s[0], nS)
        perturb_otus_P[s,remove] = 0.0
        perturb_otus[s,remove] = 0
        #ChooseRandomSubset(n, k):        
    
    temp_sums = perturb_otus_P.sum(axis=1)
    perturb_otus_P = perturb_otus_P/temp_sums[:,np.newaxis]

    return (perturb_otus,perturb_otus_P)

def main(argv):

    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="input OTU file")
    
    parser.add_argument("output_file_stub", help="output file stub")
    
    parser.add_argument('-s','--random_seed',default=23724839, type=int, 
        help=("specifies seed for numpy random number generator defaults to 23724839"))
    
    #get command line arguments  
    args = parser.parse_args()
    
    import ipdb; ipdb.set_trace()
    
    inputfile = args.input_file
    
    controlfile = args.output_file_stub+"_c.csv"
    
    treatmentfile = args.output_file_stub+"_t.csv"
    
    randomseed = args.random_seed
    
    np.random.seed(set_random_state(randomseed))
    print 'Read input OTUs from ',inputfile
    otus = p.read_csv(inputfile)

    dimO = otus.shape
    nSamples =  dimO[0]
    nOTUs = dimO[1] - 1

    #extract out numpy matrix
    otuI = otus.iloc[0:nSamples,1:nOTUs + 1]
    otu_labels = list(otuI.columns.values)    
    sample_names = otus.ix[:,0]

    input_otus = otuI.as_matrix().astype(np.float64)

    input_sums = input_otus.sum(axis=1)
    input_otus_P = input_otus/input_sums[:,np.newaxis]

    (perturb_otus,perturb_otus_P) = perturb_diversity(input_otus,input_otus_P,nSamples,nOTUs)
    
    controlSS = generate_logSS(Constants.CONTROL_A,Constants.CONTROL_B/Constants.CONTROL_A,nSamples)
    treatmentSS = generate_logSS(Constants.TREATMENT_A,Constants.TREATMENT_B/Constants.TREATMENT_A,nSamples)

    control_OTUs   = generate_MN(input_otus_P,controlSS)
    treatment_OTUs = generate_MN(perturb_otus_P,treatmentSS) 


    control_OTUs_df = generate_sample_df(control_OTUs,otu_labels,sample_names)
    control_OTUs_df.to_csv(controlfile,index=False) 

    sample_names_T = []
    for i in range(nSamples):
        sample_names_T.append(sample_names[i] + "_T")
    
    treatment_OTUs_df = generate_sample_df(treatment_OTUs,otu_labels,sample_names_T)
    treatment_OTUs_df.to_csv(treatmentfile,index=False)

    combine_df = control_OTUs_df.append(treatment_OTUs_df)
    combine_df.to_csv("Combined.csv",index=False)

    input_OTUs_df   = generate_sample_df(input_otus,otu_labels,sample_names)
    perturb_OTUs_df = generate_sample_df(perturb_otus,otu_labels,sample_names_T)
    combine_OTUs_df = input_OTUs_df.append(perturb_OTUs_df)
    combine_OTUs_df.to_csv("Combined_OTUs.csv",index=False)

    sample_C = sample_names.tolist() + sample_names_T
    perturb_C = ['False']*nSamples
    perturb_T = ['True']*nSamples 
    
    perturb = perturb_C + perturb_T 

    df = p.DataFrame(perturb, index=None, columns=['Perturb'])
    df['Samples'] = sample_C

    cols = df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df = df[cols]
    df.to_csv("Meta.csv",index=False)
    #generate treatment sample sizes

if __name__ == "__main__":
    main(sys.argv[1:])
