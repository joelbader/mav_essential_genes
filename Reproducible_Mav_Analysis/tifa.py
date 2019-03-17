#!/usr/bin/env python3
import pandas as pd 
import matplotlib.pyplot as plt
import math
import statistics as st
import statsmodels.discrete.discrete_model as smd
from statsmodels.distributions.empirical_distribution import ECDF
import statsmodels.stats.multitest as smm
import numpy as np
from scipy import stats
import time
import bisect
from scipy.special import digamma

from scipy.optimize import minimize
from scipy.special import gammaln

#import scipy.optimize
#from scipy.misc import factorial
#from scipy.spatial.distance import euclidean
#from scipy.optimize import brute
#from scipy.optimize import basinhopping

#Useful Globals Initializations
n_f = 0 #Figure index (adds one after each figure plot)
def main():
    global n_f

    ###Edit this section to adapt to a new dataset
    #Parameters
    output_csv_name = 'Attenuated_Mutants_MAC109.csv'
    out_fold = 'output/'
    
    #Build full_df from csv files of counts containing all samples
    csv = 'MAC109_raw_counts.csv'
    full_df = pd.read_csv(csv,dtype={'regulatory_class':str,'bound_moiety':str})

    #Names of samples to process for essential gene detection
    prefix = 'read_count ('
    samp_names = [prefix+'TnPool1)',prefix+'TnPool2)',prefix+'TnPool3)',prefix+'TnPool4)',prefix+'TnPool5)']

    #Define location of special site files
    special_site_files = ['LPI_sites.csv']
    
    ###This section should work without modification
    #Build c_df from full_df and samp_names
    c_df = pd.concat([full_df.iloc[:,0:6],full_df[samp_names]],sort=False,axis=1)
    c_df = c_df.rename(columns={'unique_identifier (locus_tag or record_id_start_end_strand)':'uid'})
    
    process_data(samp_names,c_df,special_site_files,out_fold,output_csv_name,plots_on = True) 

def process_data(samp_names,c_df,special_site_files,out_fold,output_csv_name,plots_on = False):
    np.random.seed(525600) #A seed is set so that the same set of genes is returned each time.
    #Pass the correctly formatted data to this function. This is the main function for identifying GD/GA/ES genes.
    global n_f

    #Rename columns for code reuse
    N = len(samp_names)
    name_ls = ['ins_Tn'+str(i+1) for i in range(N)]
    c_df = c_df.rename(columns=dict(zip(samp_names,name_ls)))

    #Contig names
    contig_ls = sorted(set(c_df['contig']))

    #Write to disk
    c_df.to_csv(out_fold+'Organized_raw_data.csv',index=False)

    #Load special site files
    ss_df_ls = [pd.read_csv(f) for f in special_site_files]

    #Iterate over contigs - each contig should be treated independently due to large difference in coverage (likely due to differences in copy number)
    pval_df_ls = [] #List to hold result dataframe
    for contig in contig_ls:
        comb_df = c_df[c_df['contig']==contig]

        #Process special site files
        site_df_ls = []
        union_set = set()
        for df in ss_df_ls:
            site_set = set(df[df['contig']==contig]['insertion_site'])
            inter_set = site_set.intersection(union_set)
            if inter_set:
                exit('Special sites are in multiple files (Should only appear in a single file):',inter_set)
            else:
                union_set = union_set.union(site_set)
            site_df_ls.append(comb_df[comb_df['insertion_site'].isin(site_set)].copy())

        #The default (not associated a motif - assumed equally permissible)
        def_df = comb_df[~comb_df['insertion_site'].isin(union_set)].copy()
        site_df_ls = [def_df] + site_df_ls

        #TA sites with the Dejesus motif
        #m1_df = comb_df[comb_df['insertion_site'].isin(motif1_df['insertion_site'])].copy()
    
        #Graph number of unique insertion sites for each sample (as you add more samples).
        if plots_on:
            n_f+=1
            cumulative_hits(comb_df,N,n_f,out_fold+'All_sites_cumulative'+'_'+contig+'.png')

            #Calculate approximate number of permissible insertion sites
            n_f+=1
            N_est1_def = cumulative_hits(site_df_ls[0],N,n_f,out_fold+'Default_sites_cumulative'+'_'+contig+'.png')
            n_f+=1
            N_est1_m1 = cumulative_hits(site_df_ls[1],N,n_f,out_fold+'Motif1_sites_cumulative'+'_'+contig+'.png')

        #Plot histogram of WT genes
        #if plots_on:
        #    n_f+=1
        #    plt.figure()
        #    data = site_df_ls[0]['ins_Tn1']
        #    weights = len(data)*[(1/len(data))]
        #    plt.hist((1/np.log(10))*np.log1p(data),bins=100,weights=weights)
        #    plt.xlabel('Log$_{10}$(Read Count + 1)')
        #    plt.ylabel('Fraction of TA sites')
        #    plt.savefig(out_fold+'Read_Distribution_Tn1'+'.png')

        #    n_f+=1
        #    plt.figure()
        #    data = site_df_ls[0][wt_ind_def]
        #    weights = len(data)*[(1/len(data))]
        #    plt.hist((1/np.log(10))*np.log1p(data),bins=100,weights=weights)
        #    plt.xlabel('Log$_{10}$(Read Count + 1)')
        #    plt.ylabel('Fraction of TA sites')
        #    plt.savefig(out_fold+'Read_Distribution_Tn1_filtered'+'.png')
        #    
        #    plt.show()
        #    exit()

        #Calculate pvalues for attenuation
        #Choose one sample as a reference - take the middle 50% of read counts. Use these names as the wildtype in the other samples.
        filter_low = 0.4
        filter_high = 0.85
        att_pval_ls = []
        ess_pval_ls = []
        for name in name_ls:
            #Use left and right ecdf to compute random smoothed attenuation pvalues for each site. Pool pvalues at each site first. Then pool over all sites in the gene second.
            print('Processing Sample: '+name)

            ref_name_ls = list(set(name_ls) - set([name])) #All samples except "name"
            #print(ref_name_ls)

            att_pval_site_ls = []
            ess_pval_site_ls = []
            for df in site_df_ls:
                #Identify mutants without a phenotype (WT)
                rank_avg_sr = df[ref_name_ls].rank(axis=0).mean(axis=1)
                wt_ind = get_wt_ind(rank_avg_sr,filter_low,filter_high)

                att_pval_sr = get_pval_sr(df[name],wt_ind)
                ess_pval_sr = essnb_pval(df[name],wt_ind,plots_on = False)

                att_pval_site_ls.append(att_pval_sr)
                ess_pval_site_ls.append(ess_pval_sr)

                #Plot rank distribution
                #n_f += 1
                #plt.figure(n_f)
                #weights = len(rank_avg_sr)*[(1/len(rank_avg_sr))]
                #plt.hist(rank_avg_sr,bins=100,weights=weights)
                #plt.savefig(out_fold+'Mean_rank_distribution'+'.png')

                #print('Motif1 Sites:')
                #att_pval_m1_sr = get_pval_sr(m1_df[name],wt_ind_m1)
                #ess_pval_m1_sr = essnb_pval(m1_df[name],wt_ind_m1,plots_on = False)
        
            temp_df = pd.concat(att_pval_site_ls).sort_index()
            temp_df.name = name + '_att'
            att_pval_ls.append(temp_df)

            temp_df = pd.concat(ess_pval_site_ls).sort_index()
            temp_df.name = name + '_ess'
            ess_pval_ls.append(temp_df)

        names_att = [name + '_att' for name in name_ls]
        names_ess = [name + '_ess' for name in name_ls]

        #Fitness Calculation
        fit_df_ls = []
        for df in site_df_ls:
            fit_df_ls.append(relative_fitness(df,N))
        fit_comb = pd.concat(fit_df_ls).sort_index()

        df = pd.concat([comb_df[['uid','product']],fit_comb] + att_pval_ls + ess_pval_ls,axis=1)

        #Remove sites which don't contain any features:
        df = df[~df['uid'].isna()]
        
        #Remove sites which contain more than one feature.
        func = lambda x: x.find(';') == -1
        df = df[df['uid'].map(func)]

        print("Pooling p-values.")
        #Compute Stouffer pooling for each TA site
        df['Z_att'] = np.sum(stats.norm.ppf(df[names_att]),axis=1)/(pow(len(names_att),0.5))
        df['Z_ess'] = np.sum(stats.norm.ppf(df[names_ess]),axis=1)/(pow(len(names_ess),0.5))

        df.to_csv(out_fold + 'Testing_'+contig+'_df.csv')
        #exit()

        #Parameters for removal of TA sites (protection against biased sites)
        minTArm = 0*1
        fracTArm = 0*0.1

        row_ls = []
        for name,g in df.groupby(['uid'],sort=False):
            L = len(g.index)
            if L > minTArm:
                #Attenuate pvalue
                Zscore = sum(g['Z_att'])/pow(L,0.5)
                if Zscore < 0:
                    #Remove bottom 10% or at least 1 site
                    num_rm = max(minTArm,round(fracTArm*L)) #Number of TA sites to remove
                    Zscore = sum(g['Z_att'].sort_values()[num_rm:L])/pow(L-num_rm,0.5)
                    p_pooled = stats.norm.sf(abs(Zscore))*2
                else:
                    p_pooled = stats.norm.sf(abs(Zscore))*2
                
                #Essential pvalue 
                num_rm = max(minTArm,round(fracTArm*L)) #Number of TA sites to remove
                p_ess = stats.norm.cdf(g['Z_ess'].sort_values()[num_rm:L])
                p_pooled_ess = tpm_pval(p_ess,t=0.5)

                #Fitness
                rel_fit = g['Relative_fitness'].median()
                #Combine
                row_ls.append([name,g['product'].iloc[0],rel_fit,p_pooled,(Zscore < 0),p_pooled_ess])
            else:
                pass
                #Combine
                #row_ls.append([name,g['product'].iloc[0],np.nan,np.nan,False,np.nan])

        pval_df = pd.DataFrame(row_ls,columns=['Gene','product','Relative Fitness','pvalue','Stouffers GD','pvalue_ess'])
        pval_df_ls.append(pval_df)

        if plots_on:
            n_f+=1
            plt.figure(n_f)
            plt.hist(np.log2(pval_df['Relative Fitness'][pval_df['Relative Fitness']>0].astype(np.float64)),bins=20, color='r')
            plt.xlabel('log2(Relative Fitness)')
            plt.ylabel('Number of genes (N='+str(len(pval_df['Relative Fitness'])) + ')')
            plt.savefig(out_fold+'Relative_fitness_hist'+'_'+contig+'.png')
            #plt.show()

    #Concat the 3 results for each contig together
    fc_thresh = 1.5 #Relative fitness threshold to be called GD or GA.
    out_df = pd.concat(pval_df_ls)
    rej,padj_ls = smm.multipletests(out_df['pvalue'], alpha=1, method='b')[0:2]
    out_df['GD'] = rej & out_df['Stouffers GD'] & (out_df['Relative Fitness'] < 1/fc_thresh) #Stouffers GD gives direction
    out_df['GA'] = rej & ~out_df['Stouffers GD'] & (out_df['Relative Fitness'] > fc_thresh)

    if (sum(out_df['GD'] & (out_df['Relative Fitness'] > 1)) > 0) | (sum(out_df['GA'] & (out_df['Relative Fitness'] < 1)) > 0):
        print('###########################################')
        print('The following genes are statistically significant but two distinct measures of the fitness are in conflict. This can happen for very long proteins where a large part of the protein is unimportant for the phenotype of the mutant.')
        if sum(out_df['GD'] & (out_df['Relative Fitness'] > 1)) > 0:
            print(out_df[out_df['GD'] & (out_df['Relative Fitness'] > 1)])
        if sum(out_df['GA'] & (out_df['Relative Fitness'] < 1)) > 0:
            print(out_df[out_df['GA'] & (out_df['Relative Fitness'] < 1)])
        print('###########################################\n')

    out_df['Adjusted p-value'] = padj_ls
    print('{0} hypotheses rejected for GD, {1} hypotheses rejected for GA, out of {2} genes'.format(sum(out_df['GD']),sum(out_df['GA']),len(rej)))

    #Essential Genes
    rej,padj_ls = smm.multipletests(out_df['pvalue_ess'], alpha=0.01, method='fdr_bh')[0:2]
    #rej,padj_ls = smm.multipletests(out_df['pvalue_ess'], alpha=1.0, method='b')[0:2]
    out_df['Essential'] = rej
    out_df['Adjusted p-value (Ess)'] = padj_ls
    print('{0} hypotheses rejected for Essential Test out of {1} genes'.format(sum(rej),len(rej)))

    print('Number of genes identified as "ES" but not "GD" (can be caused by differences in statistical cut-offs, proteins with non-essential and essential domains, and in very rare cases (very lucky/unlucky draws during p-value calculation)). By default these will be assigned "NE":',sum(out_df['Essential'] & ~out_df['GD']))
    if sum(out_df['Essential'] & ~out_df['GD']) > 0:
        print(out_df[out_df['Essential'] & ~out_df['GD']])

    #Label groups - for ease of interpretation
    out_df['Group'] = 'NE'
    out_df.loc[out_df['GD'], 'Group'] = 'GD'
    out_df.loc[out_df['GD'] & out_df['Essential'],'Group'] = 'ES' #This resolves conflicts between GD and ES (which should be rare). Conflicts will still be reported in the True/False columns.
    out_df.loc[out_df['GA'], 'Group'] = 'GA'

    #Output Dataframe
    out_df.index = out_df['Gene']
    out_df = out_df.drop(['Gene','Stouffers GD'],axis=1)
    out_df = out_df.rename({'product':'Annotation','Relative Fitness':'Estimated Fitness','pvalue':'pval (GD or GA)', 'pvalue_ess':'pval (Essential)', 'Adjusted p-value':'padj (GD or GA)','Adjusted p-value (Ess)':'padj (Essential)'},axis=1)
    out_df = out_df[['Group','GD','GA','Essential','Estimated Fitness','pval (GD or GA)','padj (GD or GA)','pval (Essential)','padj (Essential)','Annotation']] #Change order of columns
    out_df.to_csv(out_fold+output_csv_name)

def comb_pval(sr):
    _,p_pooled = stats.combine_pvalues(sr,method='stouffer')
    return p_pooled 

def tpm_pval(p_sr, t=0.5):
    p = np.array(p_sr) #Convert to numpy for speed

    L = len(p)
    lnw = sum(np.log(p[p < t]))
    #print('lnw:',lnw)
    #print('t:',t)

    #w = (p_sr[p_sr < t]).product(skipna=False) #pandas is much slower than numpy

    p_comb = 0.0
    for k in range(1,L+1):
        innersum = 0.0
        if lnw <= k*math.log(t):
            A = k*math.log(t) - lnw
            s = np.array(range(k))
            innersum = np.exp(lnw)*sum( np.exp(np.log(A)*s - gammaln(1+s)) )
        else:
            innersum = pow(t, k)

        r1 = gammaln(1+L) - gammaln(1+L-k) - gammaln(1+k)
        if 1-t < pow(10,-50):
            p_comb += math.exp(r1)*pow(1-t,L-k) * innersum
        else:
            r2 = (L-k)*math.log(1-t)
            p_comb += math.exp(r1+r2) * innersum

    return p_comb

def snb_pval(sr,t,p,r):
    pval_sr = pd.Series(np.random.uniform(t + (1-t)*stats.nbinom.cdf(sr-1,r,p),t + (1-t)*stats.nbinom.cdf(sr,r,p)))
    pval_sr.index = sr.index
    pval_sr.name = sr.name
    pval_sr[sr == 0] = np.random.uniform(0,t + (1-t)*stats.nbinom.cdf(sr[sr == 0],r,p))

    return pval_sr

def essnb_pval(sr,wt_ind,plots_on = False,scale=0.05):
    #Input: Series of read counts (sr) and wt_ind (boolean series for selection of WT mutants)
    #Optional Input: plots_on = True if you wish to produce a plot of the fit parameters against the data, scale is the relative fitness threshold for a mutant to be declared essential.
    #Output: Pandas Series giving the p-values of each insertion site.
    #Fit zero-inflated NB to WT data
    global n_f
    data = sr[wt_ind].tolist()
    ecdf1 = ECDF(data)

    t,p,r = fit_snb(data)

    #print('zero-inflated negative binomial parameters: t='+str(t)+'  mu='+str(r/p - r) + '   a='+str(1/r))
    if plots_on:
        n_f+=1
        plt.figure(n_f)
        plt.plot(sorted(data),ecdf1(sorted(data)),'b',label='ecdf')
        plt.plot(sorted(data),sorted(snb_cdf(data,(t,p,r))),'g',label='NB fit')
        plt.xlabel('Read Count')
        plt.ylabel('Fraction of mutants less than or equal to Read Count')
        plt.legend()
        plt.show()

    #Rescale NB (mean only, leaving dispersion the same)
    tscaled = t
    pscaled = r / (r + scale*(r/p - r))
    rscaled = r
    #For each element of sr, compute one-sided pvalue given NB essential distribution. Use smoothing trick as NB is discrete. 
    pval_sr = snb_pval(sr,tscaled,pscaled,rscaled)
    pval_sr.name = pval_sr.name

    return pval_sr

def get_wt_ind(sr,leftq,rightq):
    #Input: sr a pandas series containing numeric values, leftq = lower quantile, rightq = upper quantile
    #Output: Boolean pandas series.
    #Selects the indices of the values between the leftq quantile and rightq quantile
    #Use the mean of the ranks and return true for the middle fraction (between the two quantiles) rounded "inward".
    left_ind = math.ceil(leftq * len(sr))
    right_ind = math.ceil(rightq * len(sr))
    sr_temp = sr.sort_values().reset_index(drop=True)
    left_val = sr_temp[left_ind]
    right_val = sr_temp[right_ind]
    wt_ind = (sr < right_val) & (sr > left_val)
    return wt_ind

def relative_fitness(def_df,N):
    fit_df = pd.DataFrame()
    for i in range(N):
        fit_df['norm_'+str(i+1)] = def_df['ins_Tn'+str(i+1)].copy()/sum(def_df['ins_Tn'+str(i+1)])
    fit_df['Average_norm'] = fit_df.mean(axis=1)
    norm = fit_df['Average_norm'][fit_df['Average_norm'] > 0].median(axis=0) #Should represent a WT mutant
    fit_df['Relative_fitness'] = fit_df['Average_norm']/norm
    return fit_df['Relative_fitness']

def get_pval_sr(sr,wt_bool_sr):
    #Input: sr = Series of read counts, wt_bool_sr = List of indices assumed to come from WT distribution.
    #Output: Continuous smoothed pvalues
    def get_pvalue(ind,sr,wt_bool_sr,s_count_ls):
        count = sr[ind]
        if ~wt_bool_sr[ind]:
            #This mutant not used in WT list
            pval_lb = (bisect.bisect_left(s_count_ls,count)) / (len(s_count_ls)+1)
            pval_ub = (bisect.bisect_right(s_count_ls,count) + 1) / (len(s_count_ls)+1)
        else:
            #This mutant is already present in WT list
            pval_lb = (bisect.bisect_left(s_count_ls,count)) / (len(s_count_ls))
            pval_ub = (bisect.bisect_right(s_count_ls,count)) / (len(s_count_ls))
        pval = np.random.uniform(pval_lb,pval_ub)
        return(pval)

    sr_index = pd.Series(sr.index)
    s_count_ls = sorted(sr[wt_bool_sr].tolist())

    pval_sr = sr_index.apply(get_pvalue,args=(sr,wt_bool_sr,s_count_ls))
    pval_sr.index = sr.index
    pval_sr.name = sr.name

    return pval_sr

def snb_cdf(x,snbparams):
    t = snbparams[0]
    p = snbparams[1]
    r = snbparams[2]
    return(t + (1-t)*stats.nbinom.cdf(x,r,p))

def unique(seq, idfun=None): 
    # order preserving version of "set()"
   if idfun is None:
       def idfun(x): return x
   seen = {}
   result = []
   for item in seq:
       marker = idfun(item)
       if marker in seen: continue
       seen[marker] = 1
       result.append(item)
   return result

def cumulative_hits(c_df,N,n_f,filename):
    hits_ls=[None]*N#Number of hits in each sample
    cum_hits_ls=[None]*N#Number of new hits in all samples up to index (cumulative new hits)
    cum_hits = c_df['ins_Tn1'] #Initialized to first column for convenience, could have used all zeros.
    frac_ls = [None]*N
    cum_frac_ls = [None]*N

    for i in range(N):
        cum_hits = (cum_hits | c_df['ins_Tn'+str(i+1)])
        hits_ls[i] = sum(c_df['ins_Tn'+str(i+1)]>0)
        frac_ls[i] = sum(c_df['ins_Tn'+str(i+1)]>0) / len(cum_hits)
        cum_hits_ls[i] = sum(cum_hits>0)
        cum_frac_ls[i] = sum(cum_hits>0)/len(cum_hits)

    plt.figure(n_f)
    samp_ind = list(range(1,N+1))
    plt.bar(samp_ind,frac_ls,label='Per Sample')
    plt.plot(samp_ind,cum_frac_ls,'k',label='Cumulative')
    plt.xlabel('Sample number')
    plt.ylabel('Fraction of unique TA-sites (per sample or cumulative)')
    plt.ylim(0,1)
    plt.legend(loc='upper left')
    plt.savefig(filename)
    #plt.show()

    return(cum_hits_ls[-1])

#def fit_nb(r1):
#    #Fits negative binomial with Maximum-likelihood estimation via brent's method
#    r = scipy.optimize.brentq(grad_r,0.00001,1000,args=(np.array(r1),))
#    p = len(r1)*r / (len(r1)*r + np.sum(r1))
#
#    return p,r

#def grad_r(r,data):
#    #Negative of log likelihood but removed k! scaling term
#    N = len(data)
#
#    p = N*r / (N*r + np.sum(data))
#    grad_r = N*math.log(p) + np.sum(digamma(data+r)) - N*digamma(r)
#    return grad_r

def fit_snb(data):
    data_np = np.array(data)
    t0 = 1 - np.count_nonzero(data_np)/len(data_np) #t has to be less than this - use it as a starting point
    x0=np.array([t0,0.1,0.1])

    #Compute ecdf at each value of data - store as vector
    #ecdf_data = pd.Series(data).rank(method='max')/len(data)
    #bounds=[(0.00001,0.9999),(0.00001,0.9999),(0.00001,None)]
    bounds=[(0.000001,0.9999),(0.00001,0.9999),(0.00001,None)]
    res = minimize(nlogl_reduced,x0=x0,method='L-BFGS-B',args=(data_np),bounds=bounds,jac=grad_fun)
    #Basin Hopping
    #min_args = {'method':'L-BFGS-B','args':(data,ecdf_data),'bounds':bounds}
    #res = basinhopping(sq_loss,x0=x0,niter=20,T=1.0,stepsize=0.1,minimizer_kwargs=min_args)

    x = res.x.tolist()
    return x[0],x[1],x[2]

def nlogl_reduced(x,data):
    #print('eval point:',x)
    #Negative of log likelihood but removed k! scaling term
    t = abs(x[0])
    p = abs(x[1])
    r = abs(x[2])
    data_nz = data[np.nonzero(data)] #You could move this for speed?
    N = len(data)
    n = N - len(data_nz)
    sumk = np.sum(gammaln(data_nz+r)) + math.log(1-p)*np.sum(data_nz)
    #sum([gammaln(k + r) + k*math.log(1-p) for k in data_nz])
    val = -1*(n*math.log(t+(1-t)*math.pow(p,r)) + (N-n)*math.log(1-t) + sumk - (N-n)*gammaln(r) + (N-n)*r*math.log(p))
    #print('Objective:',val)
    return val

def grad_fun(x,data):
    #Negative of log likelihood but removed k! scaling term
    t = abs(x[0])
    p = abs(x[1])
    r = abs(x[2])
    data_nz = data[np.nonzero(data)] #You could move this for speed?
    N = len(data)
    n = N - len(data_nz)
    
    pr = math.pow(p,r)
    grad_t = -1*(n*(1-pr)/(t + (1-t)*pr) - (N-n)/(1-t))
    grad_p = -1*(n*(1-t)*math.pow(p,r-1)/(t + (1-t)*pr) - np.sum(data_nz)/(1-p) + (N-n)*r/p)
    grad_r = -1*(n*math.log(p)*pr*(1-t)/(t + (1-t)*pr) + np.sum(digamma(data_nz+r)) - (N-n)*digamma(r) + (N-n)*math.log(p))
    grad = np.array([grad_t,grad_p,grad_r])
    #print('Gradient:',grad)
    return grad

#Graveyard - Deprecated Functions
#def get_pval_sr(sr,wt_readcount):
#    #Input: sr = Series of read counts, wt_readcount = List of read counts assumed to be from WT
#    #Output: Continuous smoothed pvalues
#    ecdfl = ECDF(wt_readcount,side='left')
#    ecdfr = ECDF(wt_readcount,side='right')
#    def get_pvalue(count,ecdfl,ecdfr):
#        pval = np.random.uniform(ecdfl(count),ecdfr(count))
#        return(pval)
#
#    pval_sr = sr.apply(get_pvalue,args=(ecdfl,ecdfr))
#
#    return pval_sr
#
#def fit_snb_gd(data,x0=np.array([0.1,0.1,0.1])):
#    #Fits spiked negative binomial with gradient descent
#    #Deprecated - Slated for removal
#    xnew = x0
#    grad_new = grad_fun(xnew,data)
#    while True:
#        if (xnew==x0).all():
#            gamma = math.pow(10,-8)
#        else:
#            #print(xnew,xold,grad_new,grad_old)
#            gamma = np.dot((xnew - xold),(grad_new - grad_old)) / (math.pow(euclidean(grad_new, grad_old),2))
#            print('gamma:',gamma)
#        xold = xnew
#        grad_old = grad_new
#        xnew = xold - gamma*grad_new
#        grad_new = grad_fun(xnew,data)
#        if euclidean(xnew, xold) < 0.0001:
#            break
#        print(xnew)
#        #print('Objective:',nlogl_reduced(xnew,data))
#
#    xnew = xnew.tolist()
#    return abs(xnew[0]),abs(xnew[1]),abs(xnew[2])

#def fit_snb(r1):
#    #Fits spiked negative binomial with statsmodel negative binomial fitting.
#    #Deprecated - Slated for removal
#
#    #Remove the zeros
#    r1_nz = [x for x in r1 if x>0]
#    #r1_nz = r1
#    t = (len(r1) - len(r1_nz))/len(r1)
#    #Fit to remaining values
#    result = smd.NegativeBinomial(r1_nz,np.ones(len(r1_nz))).fit(method = 'bfgs')
#    mu = math.exp(result.params[0])
#    a = result.params[1]
#
#    r = 1/a
#    p = r / (mu + r)
#    print(mu,a)
#    print(p,r)
#
#    return t,r,p
#
#
#
#def fit_snb_lsq(data,x0=np.array([0.1,0.1,0.1])):
#    #Compute ecdf at each value of data - store as vector
#    ecdf_data = pd.Series(data).rank(method='max')/len(data)
#    #bounds=[(0.00001,0.9999),(0.00001,0.9999),(0.00001,None)]
#    bounds=[(0,0),(0.00001,0.9999),(0.00001,None)]
#    res = minimize(sq_loss,x0=x0,method='L-BFGS-B',args=(data,ecdf_data),bounds=bounds)
#    #Basin Hopping
#    #min_args = {'method':'L-BFGS-B','args':(data,ecdf_data),'bounds':bounds}
#    #res = basinhopping(sq_loss,x0=x0,niter=20,T=1.0,stepsize=0.1,minimizer_kwargs=min_args)
#
#    x = res.x.tolist()
#    return x[0],x[1],x[2]

#def sq_loss(snbparams,data,ecdf_data):
#    #Input: data - data points, ecdf_data - empirical cdf evaluated at data, snbparams - tuple (t,p,r) of negative binomial params.
#    #Output: Sum of squares between ecdf_data and fit spiked negative binomial cdf at each data.
#    loss = np.sum((snb_cdf(data,snbparams) - ecdf_data)**2)
#    return loss
#
#def fit_snb_custom(data,x0=np.array([0.1,0.1,0.1]),read_thresh=0):
#    #Fits spiked negative binomial with scipy minimize().
#    data_nz = [x for x in data if x>read_thresh] #Threshold for removing strongly attenuated genes - necessary to avoid underestimation of p-values.
#    data = np.array([0]*(len(data) - len(data_nz)) + data_nz)
#    #bounds=[(0.00001,0.9999),(0.00001,0.9999),(0.00001,None)]
#    bounds=[(0,0),(0.00001,0.9999),(0.00001,None)]
#    res = minimize(nlogl_reduced,x0=x0,method='L-BFGS-B',jac=grad_fun,args=data,bounds=bounds)
#    options = {'maxfev':100,'rhobeg':0.01}
#    #res = minimize(nlogl_reduced,x0=x0,method='SLSQP',jac=grad_fun,options=options,args=data,bounds=bounds)
#    #Basin Hopping
#    min_args = {'method':'L-BFGS-B','args':data,'jac':grad_fun,'bounds':bounds}
#    #res = basinhopping(nlogl_reduced,x0=x0,niter=20,T=1.0,stepsize=0.1,minimizer_kwargs=min_args)
#    x = res.x.tolist()
#
#    #Brute Force Method
#    #x = brute(nlogl_reduced,ranges=((0.001,0.2),(0.01,0.2),(0.01,2.0)),args=(data,),Ns=10,finish=None)
#    #x = x.tolist()
#
#    #print(x)
#
#    return x[0],x[1],x[2]


if __name__ == '__main__':
    main()
