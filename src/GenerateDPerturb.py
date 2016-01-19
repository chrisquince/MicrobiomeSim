import sys, getopt
import pandas as p
import numpy as np
import random as rnd
from itertools import compress
import argparse

#sample without replacement
def sample_without_replacement(freq_matrix,sizes):

    X = np.zeros(freq_matrix.shape,dtype=np.int)

    nsamples = freq_matrix.shape[0]
    nOTUs = freq_matrix.shape[1]
    for i in range(nsamples):
        elements = np.repeat(range(nOTUs), freq_matrix[i,:].astype(int))
        sample_e = np.random.choice(elements, size=sizes[i], replace=False)
        unique, counts = np.unique(sample_e, return_counts=True)
        nC = counts.shape[0]
        for j in range(nC):
            X[i,unique[j]] += counts[j]

    return X

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

def perturb_diversity(input_otus,input_otus_P,nSamples,nOTUs,dReduce):

    perturb_otus_P = np.copy(input_otus_P)
    perturb_otus   = np.copy(input_otus)
    for s in range(nSamples):
        #get list of nonzero indices
        nonzero_s = np.nonzero(input_otus_P[s,:])
        nP = len(nonzero_s[0])
        nS = int(dReduce*nP)
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

    
    #get command line arguments  
    args = parser.parse_args()
    
    #import ipdb; ipdb.set_trace()
    
    inputfile = args.input_file
    
    controlfile = args.output_file_stub+"_c.csv"
    controlfileR = args.output_file_stub+"_c_R0.csv"
    treatmentfile = args.output_file_stub+"_t.csv"
    treatmentfileR = args.output_file_stub+"_t_R0.csv"
    combinedfile = args.output_file_stub+"_ct.csv"
    combinedotusfile = args.output_file_stub+"_ct.csv"
    metafile = args.output_file_stub+"_m.csv"
    rarefile = args.output_file_stub+"_r.csv"
        
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

    rare_sizes = np.repeat(args.base_level,nSamples)

    input_rare = np.copy(input_otus)
#sample_without_replacement(input_otus,rare_sizes)
    
    #remove empty columns
    selectC = input_rare.sum(axis=0) > 0
    otu_labels_rare = list(compress(otu_labels, selectC)) 
    input_rare = input_rare[:, selectC]
    
    input_rare_df = generate_sample_df(input_rare,otu_labels_rare,sample_names)
    input_rare_df.to_csv(rarefile,index=False)

    input_sums = (input_rare.sum(axis=1)).astype(float)

    input_rare_P = input_rare.astype(float)/input_sums[:,np.newaxis]

    nROTUs = input_rare.shape[1]
    (perturb_otus,perturb_otus_P) = perturb_diversity(input_rare,input_rare_P,nSamples,nROTUs,args.diversity_reduction)
    
    controlSS = generate_logSS(args.control_a,args.control_b/args.control_a,nSamples)
    treatmentSS = generate_logSS(args.treatment_a,args.treatment_b/args.treatment_a,nSamples)

    control_OTUs   = generate_MN(input_rare_P,controlSS)
    treatment_OTUs = generate_MN(perturb_otus_P,treatmentSS) 

    control_OTUs_df = generate_sample_df(control_OTUs,otu_labels_rare,sample_names)
    control_OTUs_df.to_csv(controlfile,index=False) 

    selectCR = control_OTUs.sum(axis=0) > 0
    otu_labels_CR = list(compress(otu_labels_rare, selectCR)) 
    control_OTUsR = control_OTUs[:, selectCR]
    control_OTUs_dfR = generate_sample_df(control_OTUsR,otu_labels_CR,sample_names)
    (control_OTUs_dfR.transpose()).to_csv(controlfileR,header=False) 

    sample_names_T = []
    for i in range(nSamples):
        sample_names_T.append(sample_names[i] + "_T")
    
    treatment_OTUs_df = generate_sample_df(treatment_OTUs,otu_labels_rare,sample_names_T)
    treatment_OTUs_df.to_csv(treatmentfile,index=False)

    selectTR = treatment_OTUs.sum(axis=0) > 0
    otu_labels_TR = list(compress(otu_labels_rare, selectTR)) 
    treatment_OTUsR = treatment_OTUs[:, selectTR]
    treatment_OTUs_dfR = generate_sample_df(treatment_OTUsR,otu_labels_TR,sample_names) 
    (treatment_OTUs_dfR.transpose()).to_csv(treatmentfileR,header=False)

    combine_df = control_OTUs_df.append(treatment_OTUs_df)
    combine_df.to_csv(combinedfile,index=False)

    input_OTUs_df   = generate_sample_df(input_otus,otu_labels,sample_names)
    perturb_OTUs_df = generate_sample_df(perturb_otus,otu_labels_rare,sample_names_T)
    combine_OTUs_df = input_OTUs_df.append(perturb_OTUs_df)
    combine_OTUs_df.to_csv(combinedotusfile,index=False)

    sample_C = sample_names.tolist() + sample_names_T
    perturb_C = ['False']*nSamples
    perturb_T = ['True']*nSamples 
    
    perturb = perturb_C + perturb_T 

    df = p.DataFrame(perturb, index=None, columns=['Perturb'])
    df['Samples'] = sample_C

    cols = df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df = df[cols]
    df.to_csv(metafile,index=False)
    #generate treatment sample sizes

if __name__ == "__main__":
    main(sys.argv[1:])
