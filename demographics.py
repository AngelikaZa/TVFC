# Demographics

## Author: A Zarkali
## Last edited: 29/10/2020
## Aim: Useful functions to make group comparisons for demographics easy

#Import libraries
import pandas as pd
import numpy as np
import scipy.stats as stats
import statsmodels.api as sm
from scipy.stats import shapiro
from statsmodels.formula.api import ols
import scikit_posthocs as sp
from pandas.api.types import is_numeric_dtype


# Function to check normality for all variables across a dataframe and return a list of normally distributed variables and one of non-normally distributed ones
### INPUTS:
###### dataframe: the pandas dataframe
### OUTPUTS:
###### returns two lists: of normally distributed and non normally distributed variables in the dataframe
def normalityCheck(df): 
    ### Declare empty variables to hold column names
    NormallyDistributed = []
    NonNormallyDistributed = []
    ### Loop through all columns
    for col in df.columns:
        if is_numeric_dtype(df[col]) == True: ## Numeric check
            data = df[np.isfinite(df[col])] ## Drop NAs (the shapiro will not calculate statistic if NAs present)
            r, p = stats.shapiro(data[col]) ### If less than 0.05 non normally distributed
            if p < 0.05: 
                NonNormallyDistributed.append(col)
            else:
                NormallyDistributed.append(col)
    return(NormallyDistributed, NonNormallyDistributed)

# Function to compare all numerical variables of a dataframe across 2 or more groups:
### INPUTS: 
###### variables: the list of variables that you want to compare between groups (they will be columns of the dataframe)
###### group: the name of the variable describing the groups you want to compare as a string (will be a column of the dataframe)
###### dataframe: the pandas dataframe 
###### number_groups: integer, the number of groups you want to compare 
### OUTPUTS: 
###### returns a new dataframe with the characteristis, statistics (test and p value) and the name of the statistical test that was used
def groupCompare(variables, group, dataframe, number_groups): 
    ### Declare empty variables to hold column names
    NormallyDistributed = []
    NonNormallyDistributed = []
    statistic = []
    p_value = []
    types = []
    ### Loop through all columns of a dataframe and check normality
    for col in dataframe.columns:
        if is_numeric_dtype(dataframe[col]) == True: ## Numeric check
            data = dataframe[np.isfinite(dataframe[col])] ## Drop NAs (the shapiro will not calculate statistic if NAs present)
            r, p = stats.shapiro(data[col]) ### If less than 0.05 non normally distributed
            if p < 0.05: 
                NonNormallyDistributed.append(col)
            else:
                NormallyDistributed.append(col)
    for var in variables:
        if number_groups > 2: 
            if var in NormallyDistributed: ## Normally distributed then do ANOVA
                data=dataframe[np.isfinite(dataframe[var])]
                variable = data[var].dropna()
                comp = data[group] ### comparison of interest
                anova = ols("variable ~ C(comp)", data=data).fit() ### run anova  
                r = anova.rsquared_adj ## extract overall model adjusted r statistic
                p = anova.f_pvalue ## extract overall model p-value
                statistic.append(r)
                p_value.append(p)
                types.append("ANOVA")
            elif var in NonNormallyDistributed: ### Non normally distributed then do Kruskal Wallis
                data = dataframe[np.isfinite(dataframe[var])] 
                ### declare the three series
                v1 = data[data[group] == 0][var] 
                v2 = data[data[group] == 1][var]
                v3 = data[data[group] == 2][var]
                r,p = stats.kruskal(v1, v2, v3) ### run Kruskal wallis
                statistic.append(r)
                p_value.append(p)
                types.append("Kruskal-Wallis")
            else: ### In case any variables were labelled incorrectly
                statistic.append("NA")
                p_value.append("NA")
                types.append("NA")
        elif number_groups == 2:
            if var in NormallyDistributed: ## Normally distributed then do ttest
                data=dataframe[np.isfinite(dataframe[var])]
                v1 = data[data[group] == 1][var]
                v2 = data[data[group] == 2][var]
                r, p = stats.ttest_ind(v1, v2)
                statistic.append(r)
                p_value.append(p)
                types.append("t-test")
            elif var in NonNormallyDistributed: ### Non normally distributed then do Mann-Whitney
                data = dataframe[np.isfinite(dataframe[var])] 
                v1 = data[data[group] == 1][var]
                v2 = data[data[group] == 2][var]
                r,p = stats.mannwhitneyu(v1, v2) ### run Kruskal wallis
                statistic.append(r)
                p_value.append(p)
                types.append("Mann-Whitney")
            else: ### In case any variables were labelled incorrectly
                statistic.append("NA")
                p_value.append("NA")
                types.append("NA")    
    ### Combine results on dataframe
    results = pd.DataFrame(data=np.zeros((len(variables), 0))) # empty dataframe
    results["Variable"] = variables # variable names
    results["Statistic"] = statistic # statistic
    results["Pvalue"] = p_value # p_value
    results["Type"] = types # type of statistical test used
    return(results)

# Function to compare categorical variables between groups
### INPUTS: 
###### variables: the list of variables that you want to compare between groups (they will be columns of the dataframe)
###### group: the name of the variable describing the groups you want to compare as a string (will be a column of the dataframe)
###### dataframe: the pandas dataframe 
### OUTPUTS: 
###### returns a new dataframe with the characteristis, statistics (test and p value)
def categoricalCompare(variables, group, df):
    ## Empty variables to hold results
    statistic = []
    p_value = []
    for var in variables:
        data=df[np.isfinite(df[var])]
        tab = pd.crosstab(data[group],data[var])
        r, p, dof, exp = stats.chi2_contingency(tab)
        statistic.append(r)
        p_value.append(p)
    ### Combine results on dataframe
    results = pd.DataFrame(data=np.zeros((len(variables), 0))) # empty dataframe
    results["Variable"] = variables # variable names
    results["ChiSquare"] = statistic # statistic
    results["Pvalue"] = p_value # p_value
    return(results)

# Function to return group statistics, mean and std for a dataframe
### IMPUTS: 
###### df: the pandas dataframe
###### group: the name of the variable defining your groups (a string)
### OUTPUTS: 
###### a dataframe with mean and sd for each variable per group
def groupDemographics(df, group):
    ## Save all group mean and std for each column
    ### Group by PD_VHAny: 0 controls, 1 PD non VH, 2 PD VH
    dfMean = df.groupby(group).mean()
    dfMean["Type"] = "Mean"
    dfSTD = df.groupby(group).std()
    dfSTD["Type"] = "STD"
    result = pd.concat([dfMean,dfSTD], axis=0)
    return result