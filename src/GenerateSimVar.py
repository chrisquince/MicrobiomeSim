import sys, getopt
import pandas as p
import numpy as np

class Constants(object):
    #LOG_SCALE_POP = np.log(1.0e7)
    LOG_SCALE_POP = 0.0
    SAMPLE_SIZE   = 2000
    N_CONTROL     = 50
    N_TREATMENT   = 50
    N_CHANGE      = 25
    THRESHOLD     = -8.75
    MAX_CHANGE    = np.log(10.0)
    MIN_CHANGE    = -np.log(10.0)

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

#generate log normal relative abundances in control samples
def generate_LN_rel_abundance(means_otus_lP,std_otus_lP,nSamples,nOTUs):
    Z = np.random.normal(size=(nSamples,nOTUs))

    ZT = Z.transpose()
    
    ZU = ZT*std_otus_lP[:,np.newaxis]

    ZV = ZU.transpose()
    
    ZV = np.exp(ZV + means_otus_lP)

    sample_sizesZV = ZV.sum(axis=1)
    ZVP = ZV / sample_sizesZV[:, np.newaxis]
    return ZVP

#generate multinomial samples from matrix of relative abundances
def generate_MN(rel_abund_matrix,n_Multi_Samples):

    ra_dim = rel_abund_matrix.shape
    nSamples = ra_dim[0]
    nOTUs = ra_dim[1]

    otus = np.zeros(shape=(nSamples,nOTUs),dtype=np.int)
    for i in range(nSamples):
        otus[i,:] = np.random.multinomial(n_Multi_Samples,rel_abund_matrix[i,:])

    return otus

def generate_sample_df(control_otus,treatment_otus,cols,sample_labels):
    sampled_otus = np.vstack((control_otus,treatment_otus))

    sampled_otus_df = p.DataFrame(sampled_otus,columns=cols)

    sampled_otus_df['Id'] = sample_labels

    cols = sampled_otus_df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    sampled_otus_df = sampled_otus_df[cols]
    return sampled_otus_df


def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:s:",["ifile=","ofile=","seed="])
    except getopt.GetoptError:
        print 'test.py -i <inputfile> -o <outputfile>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'test.py -i <inputfile> -o <outputfile>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-s", "--seed"):
            seed = arg
    
    np.random.seed(set_random_state(seed))

    print 'Read input OTUs from ',inputfile
    otus = p.read_csv(inputfile)
    print otus
    
    dimO = otus.shape   
    nSamples =  dimO[0]
    nOTUs = dimO[1] - 1 
    
    #extract out numpy matrix
    otuI = otus.iloc[0:nSamples,1:nOTUs + 1]
    otus_M = otuI.as_matrix()
    
    #add pseudocount
    otus_M+=1

    #normalise rows
    otus_M = otus_M.astype(np.float64)
    sample_sizes = otus_M.sum(axis=1)
    otus_P = otus_M / sample_sizes[:, np.newaxis]
    print otus_P

    otus_logP = np.log(otus_P)
    print otus_logP

    mean_otus_lP = np.mean(otus_logP, axis=0)
    mean_otus_lPO = np.copy(mean_otus_lP)
    mean_otus_lP += Constants.LOG_SCALE_POP
 
    std_otus_lP = np.std(otus_logP,axis=0)

    print mean_otus_lP

    #generate relative abundances in controls
    controlP = generate_LN_rel_abundance(mean_otus_lP,std_otus_lP,Constants.N_CONTROL,nOTUs)

    print "\n"
    print controlP
    
    control_otus = generate_MN(controlP,Constants.SAMPLE_SIZE)
    print control_otus

    threshold = []
    for i in range(nOTUs):
        if(mean_otus_lPO[i] > Constants.THRESHOLD):
            threshold.append(i)

    print threshold

    perturb = np.random.choice(threshold, Constants.N_CHANGE, replace=False)
    print perturb

    mean_otus_lP_t = np.copy(mean_otus_lP)
    perturbO = ["False"]*nOTUs
    for i in perturb:
        randP = np.random.uniform(Constants.MIN_CHANGE, Constants.MAX_CHANGE)
        print str(randP) + "\n";
        mean_otus_lP_t[i] += randP
        perturbO[i] = "True"

    treatmentP = generate_LN_rel_abundance(mean_otus_lP_t,std_otus_lP,Constants.N_TREATMENT,nOTUs) 
    print treatmentP
    print "\n"
    treatment_otus = generate_MN(treatmentP,Constants.SAMPLE_SIZE)
    print treatment_otus
    print "\n"
    
    cols = otus.columns.tolist()
    cols.pop(0)
    sample_labels = []
    totals = Constants.N_CONTROL + Constants.N_TREATMENT
    for i in range(totals):
        sample_labels.append("Sample"+str(i))

    sampled_otus_df = generate_sample_df(control_otus,treatment_otus,cols,sample_labels)
    sampled_otus_df.to_csv(outputfile+".csv",index=False)

    meta = ["Control"]*Constants.N_CONTROL + ["Treatment"]*Constants.N_TREATMENT 

    meta_df = p.DataFrame({'Id': sample_labels, 'Group': meta})

    meta_df.to_csv(outputfile+"_meta.csv",index=False,cols=["Id","Group"],engine='python')
    
    otus_mean_df = p.DataFrame({'OTU': cols, 'Original': mean_otus_lPO, 'Control': mean_otus_lP,'Treatment':mean_otus_lP_t,'Std':std_otus_lP,'Perturb':perturbO})

    otus_mean_df.to_csv(outputfile+"_abund.csv",index=False,cols=["OTU","Original","Control","Treatment","Std","Perturb"],engine='python')

if __name__ == "__main__":
    main(sys.argv[1:])
