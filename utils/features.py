from numpy import array,exp,amin,amax,less,linspace,argmin
from sklearn import mixture
from scipy.signal import argrelextrema

################################################## estimate cutoffs between bimodal populations
def estimate_cutoff(channel) :
    X = channel.reshape(-1, 1)
        
    gmm = mixture.BayesianGaussianMixture(
        n_components=3, covariance_type='spherical',

        weight_concentration_prior=1e-2,
        mean_precision_prior=1e-2,

        covariance_prior= array([5.0]),
        mean_prior= array([0.0]),

        max_iter=5000,
        tol=2e-2

    ).fit(X)
    
    x = linspace(amin(channel),amax(channel),100)
    px = exp(gmm.score_samples(x.reshape(-1, 1)))
    
    min_index, = argrelextrema(px,less)
    xminima = x[min_index]
    pminima = px[min_index]
    
    if len(xminima)>0 :
        cutoff = xminima[argmin(pminima)]

    else : # fwhm
        min_index, = argrelextrema(abs(px-amax(px)/4),less)
        left,right = x[min_index]
        cutoff = right
        
        if sum(px[x<cutoff]*x[x<cutoff])/sum(px[x<cutoff]) > 3.0 :
            cutoff = left

    return cutoff