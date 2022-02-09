import pandas as pd
import numpy as np

def adjustsemantic(s, f):
    return (2-s-abs(s))/2 + (s+abs(s)-1)*f

def takechrbase(ref):
    size = ref.shape[0]
    chrbase = np.empty(size, dtype='object')
    for j in np.arange(0,size,1):
        chrbase[j] = "chr"+str(ref.chrom[j])+"."+str(ref.base[j])
    return chrbase

def comparechrbase(ref, input_file):
    k = 0
    size = ref.shape[0]
    freqc = np.zeros(size)
    for i in ref.chrBase:
        for j in np.arange(0,input_file.shape[0],1):
            if i == input_file.chrBase[j]:
                freqc[k] = input_file.freqC[j]
                k += 1
    return freqc

def computer(input_file, output_file):
    # ---------- Reading data: MethylDackel output in tabular format ------------------    
    #fname = fd.askopenfilename()
    input_file = pd.read_csv(input_file, header = 'infer', sep = "\t")
    # ---------- Building the reference table with chromosomes, base positions, weights, etc.
    chrom = [5,8,19,22,17,20,2,2,3,11,11,2,2]
    base = [1279643,9903378,15232074,19723434,31559925,58817319,80304673,97724391,122577739,128694264,131911273,181458175,240820204]
    gene = ["TERTmin","LINC0059_1","EPHX3","GP1BB1minus","MIR193_12","MIR296","LRRTM3","ZAP70_16","PARP15_2","FLI11","NTM14","ITGA4_3","KIF1A22"]
    strand = ["N","N","N","P","N","N","N","P","P","P","P","P","N"]
    sm = [-1, 1 , 1, -1, 1, -1, 1, 1, 1, 1, 1, 1, 1]
    weight = [16.6,-1.034,-0.034,3.896,0.206,7.736,1.579,1.778,1.526,-3.165,-0.085,1.716,1.35]
    K = [-4.304,np.NaN,np.NaN,np.NaN,np.NaN,np.NaN,np.NaN,np.NaN,np.NaN,np.NaN,np.NaN,np.NaN,np.NaN]
    #chrBase = np.empty(13,str)
    column_names = ["gene","chrom","base","strand","sm","weight","K"]#,"chrBase"]
    d = {'col1':gene, 'col2':chrom,'col3':base,'col4':strand,'col5':sm,'col6':weight,'col7':K}#, 'col8':chrBase}
    ref = pd.DataFrame(d)
    ref.columns=column_names
    ref.index=gene
    # ----------- Table lookup seaching matches for the base positions from the reference for each gene --------
    #ref.loc[:, 'chrBase'] = takechrbase(ref)
    ref = ref.assign(chrBase = takechrbase(ref))
    ref = ref.assign(freqC = comparechrbase(ref, input_file))
    ref.freqC.fillna(0)
    ref = ref.assign(freq100 = ref.freqC/100)
    ref = ref.assign(freqadj = adjustsemantic(ref.sm,ref.freq100))
    result = np.inner(ref.weight, ref.freqadj) + ref.K[0]

    # ----------- write results ----------#
    file = open(output_file,"w") 
    #file.write(display(ref))
    with open(output_file, 'a') as f:
        dfAsString = ref.to_string(header=True, index=False)
        f.write(dfAsString)
    file.close()
    file = open(output_file,"a") 
    file.write("\n\nScore for patient-file "+ str(fname)+ ": "+ str(result))
    file.close()
    
import argparse
import os

def existentFile(string):
    if os.path.isfile(string):
        return string
    else:
        msg = "could not find file %s" % string
        raise argparse.ArgumentTypeError(msg)

parser = argparse.ArgumentParser()

parser.add_argument('--inputFile', action='store', dest='inputFile', required=True, type=existentFile,
                    help='input methyldackel csv file')
parser.add_argument('--outputFile', action='store', dest='outputFile', required=True, type=existentFile,
                    help='output txt file')

args = parser.parse_args()

inputFile = args.inputFile
outputFile = args.outputFile

computer(inputFile, outputFile)
