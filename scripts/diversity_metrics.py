import sys,csv
import pandas as pd
import numpy as np


def hsn(df):
    x = np.array(df.freq)
    if len(x) > 0:
        return(sum(-1*x*np.log(x))/len(x))
    else:
        return(0)

total_sites = 2600
def site_hsn(site):
    """
    Calculate the entropy at any specific site

    """
    #create a dict for mutants
    mut_dict = dict()
    #add the non reference sequences
    for alt in site.alt:
        mut_dict[alt] = float(site[site['alt'] == alt].freq)

    #add the reference information
    mut_dict[list(site.ref)[0]] = 1 - sum(list(mut_dict.values()))

    #compute entropy using the values
    x = np.array(list(mut_dict.values()))
    return(-1*sum(x*np.log(x)))

def mean_hsn(df,total_sites=2600):
    """
    Calculates the average entropy across a list of sites

    """
    hsn_data = []
    for site in list(set(df.pos)):
        hsn_data.append(site_hsn(df[df['pos'] == site]))
    return(sum(hsn_data)/total_sites)

def L1_norm(df1,df2):

    #common snvs
    p = np.array(df1[df1.pos_alt.isin(df2.pos_alt)].freq)
    q = np.array(df2[df2.pos_alt.isin(df1.pos_alt)].freq)
    L1 = sum(abs(p - q))

    #discordant sites
    p = np.array(df1[~df1.pos_alt.isin(df2.pos_alt)].freq)
    L1 += sum(p)
    q = np.array(df2[~df2.pos_alt.isin(df1.pos_alt)].freq)
    L1 += sum(q)

    return(L1)

col_names = ['sample','hsn']

#create matrix for sample distance
number_of_samples = len(snakemake.input)
matrix = np.zeros(shape=(number_of_samples,number_of_samples))

# create lists for distance data
samples = []
sample_names = []

#write out sample specific diversity data
with open(snakemake.output[0],'w') as out_file:
    writer = csv.writer(out_file)
    writer.writerow((col_names))

    for snakefile in snakemake.input:

        with open(snakefile) as file_name:
            file = pd.read_csv(file_name)

            #add df to sample list
            samples.append(file)
            sample_names.append(file_name.name.split("/")[2].split("_")[0])

            #write diversity data
            row_data = [file_name.name.split("/")[2].split("_")[0],mean_hsn(file)]
            writer.writerow(row_data)


#write out comparative diversity data_path

for i in range(0,number_of_samples):
    for j in range(0,number_of_samples):
        matrix[i,j] = L1_norm(samples[i],samples[j])

np.savetxt(snakemake.output[1], matrix, delimiter=",",header=','.join(sample_names))
