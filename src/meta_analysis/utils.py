import pandas as pd


def compute_meta_data(dataf):
    """
    This function will compute the standardization computations.
    Overall, see Hedges & Olkin, 1985 for reference

    :param dataf: source data frame
    :return: dataframe with added columns per computation
    """
    return (dataf
            .assign(vi_SE2           = lambda d:  1/d['weight factor'])
            # compte unbiase z value in case of low sample soe, See Nakagawa & Cuthill, 2007
            .assign(standardized_ESU = lambda d:  d['standardized es'] * (1-(3/(4*(d['n'])-9))))
            .assign(wES              = lambda d:  d['standardized es'] * d['weight factor'])
            .assign(wES2             = lambda d:  d['weight factor']   * d['standardized es'] * d['standardized es'])
            .assign(w2               = lambda d:  d['weight factor']   * d['weight factor']))


def compute_random_effects(dataf, levels=['specie_str','structure_renamed']):
    """
    This function will compute the random effects per level of interest
    :param dataf: source dataframe
    :return: dataframe groupedby levels of interest with added computations
    """
    return (dataf
            .groupby(levels)
            .apply(lambda d: pd.Series({
            "sum_wi":          d['weight factor'].sum(),
            "sum_wi2":         d['w2'].sum(),
            "k":               d['w2'].count(), #num studies
            "wxES":            d['wES'].sum(),
            "wxES2":           d['wES2'].sum()
            })))


def compute_random_variance(dataf):
    """
    Compute random variance, according to random effect model.
    he Q value is a measure of the dispersion of the effect sizes.
    This measure follows the chi square distribution with k-1 degrees of freedom,
    where k is the total number of effect sizes.

    v0 is the variance due to intrinsic sampling error, according to random effect models

    See Lipsey & Wilson, 2001; Nakagawa & Cuthill, 2007

    :param dataf: source dataframe
    :return: return Q and v0 computations
    """
    return (dataf
            .assign(Q  = lambda d:  abs(d['wxES2']-(d['wxES']**2/d['sum_wi'])))
            # For struture with only 1 effect size, the denominator of v0 will automatically
            # be equal to 0. hence v0 will not be computed.
            .assign(v0 = lambda d:  (d['Q']-(d['k']-1))/(d['sum_wi']-
                                                        (d['sum_wi2']/d['sum_wi'])))
           )

def zero_if_negative(dataf):
    """
    Corrects v0 if negative value

    Reference missing!

    :param dataf: source dataframe
    :return: corrected dataframe
    """
    num = dataf['v0']
    num[num < 0] = 0
    dataf['v0'] = num
    return dataf

def apply_corrections(dataf):
    """
    This function applies the correction required by a random effect model.

    :param dataf: source dataframe
    :return: corrected dataframe
    """
    return (dataf
            .assign(v0_plus_vi       = lambda d:  d['vi_SE2']+d['v0'])
            .assign(wi_corr          = lambda d:  1/d['v0_plus_vi'])
            .assign(wxES_corr        = lambda d:  d['wi_corr']*d['standardized_ESU'])
            .assign(wxESsq_corr      = lambda d:  d['wi_corr']*(d['standardized_ESU']**2))
            .assign(SE_corr          = lambda d:  (1/(d['wi_corr'])**0.5))
            .fillna(0) # NaN are replaced by 0
            .replace(np.inf, 0)
)

def calculate_constants(dataf):
    """
    Compute constants

    Missing reference!

    :param dataf:
    :return: Corrected dataframe
    """
    return (dataf
            .groupby(['specie_str', 'structure_renamed'])
            .apply(lambda d: pd.Series({
                        "sum_wxES_corr":        d['wxES_corr'].sum(),
                        "sum_wxESsq_corr":      d['wxESsq_corr'].sum(), #num studies
                        "sum_wi_corr":          abs(d['wi_corr']).sum(),
             }))
             .fillna(0) # NaN are replaced by 0
             .replace(np.inf, 0)
           )


def I2(Q,dfg):
    """
    Computes the I^2 value, ie, percent of variation due to heterogeneity rather than chance
    By convention, Q = 0 if Q < k-1, so that the precision of a random effects summary estimate
    will not exceed the precision of a fixed effect summary estimate
    See Higgins & Thompson 2002; DOI: 10.1002/sim.1186
    :param Q: Q statistics (measure of the dispersion of the effect sizes)
    :param dfg: degress of freedom
    :return: I2 statistic
    """
#     if pd.isnull((Q < dfg)):
#         I2 = 0
#         return I2
#     else:
    I2=((Q-dfg)/Q)*100
    return I2

def calculate_mean_se(dataf):
    """
    Calculate mean and standard error of effect sizes

    For each variable and its different levels, one can calculated the mean effect size, 9
    5% confidence intervals (CI) and zscore value using eqs. 16, 17, 18 and 19. S

    See Nakagawa & Cuthill, 2007 for more information on this topic.

    Note that there are other proposed methods to compare these values.

    :param dataf: source dataframe
    :return: dataframe with added computations
    """
    return (dataf
            .assign(ES_mean       = lambda d: d['sum_wxES_corr']/d['sum_wi_corr'])
            .assign(SE_mean       = lambda d: (1/d['sum_wi_corr'])**0.5)
            .assign(z             = lambda d: d['ES_mean']/d['SE_mean'])
            .assign(high_CI       = lambda d: d['ES_mean']+(1.96*d['SE_mean']))
            .assign(low_CI        = lambda d: d['ES_mean']-(1.96*d['SE_mean']))
           # .assign(k_val         = lambda d: d['w2'].count())#num studies
            .assign(I_val         = lambda d: I2(d['Q'],d['k']-1))
           )

def calculate_number_es(dataf):
    """
    Return count of effect size per level of interest
    :param dataf: source dataframe
    :return: corrected dataframe
    """
    return (dataf
            .groupby(['specie_str', 'structure_renamed'])
            .apply(lambda d: pd.Series({
                        "k_val": d['w2'].count(),
             }))
           )