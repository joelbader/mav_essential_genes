#!/usr/bin/env python3

from tifa import process_data
import pandas as pd

#Useful Globals Initializations
n_f = 0 #Figure index (adds one after each figure plot)
def main():
    global n_f

    #Parameters
    output_csv_name = 'Attenuated_Mutants_H37Rv.csv'
    out_fold = 'output/'
    
    #Read in csv files of counts
    csv = 'H37Rv_raw_counts.csv'
    
    #Initializations
    c_df = pd.read_csv(csv,dtype={'regulatory_class':str,'bound_moiety':str})

    samp_names = c_df.columns.values[4:]

    #Rename columns for ease of use and compatibility
    c_df = c_df.rename(columns={'unique_identifier (locus_tag or record_id_start_end_strand)':'uid'})

    special_site_files = ['LPI_sites.csv']

    process_data(samp_names,c_df,special_site_files,out_fold,output_csv_name,plots_on = False) 

if __name__=='__main__':
    main()
