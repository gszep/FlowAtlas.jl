from numpy import gradient,pi,argsort,array,argmax,sqrt,sort,exp,amin,amax,less,greater,linspace,argmin,quantile,histogram,isclose
from sklearn import mixture
from scipy.signal import argrelextrema

from matplotlib.pyplot import plot,figure,show,axvline,fill_between

################################################## estimate cutoffs between bimodal populations
def estimate_cutoff(channel,cutoff,npoints=200,ntile=0.95) :
    
    # fit gaussian mixture mpodel
    mixture_model = mixture.BayesianGaussianMixture(
        n_components=3, covariance_type='spherical',

        weight_concentration_prior=1e-5,
        mean_precision_prior=1e-2,

        covariance_prior= array([5.0]),
        mean_prior= array([0.0]),

        max_iter=5000,
        tol=1e-5

    ).fit(channel.reshape(-1, 1))

    # sample probability density
    x = linspace(amin(channel),amax(channel),npoints)
    px = exp(mixture_model.score_samples(x.reshape(-1, 1)))

    ws = mixture_model.weights_.flatten()
    mus = mixture_model.means_.flatten()
    sigmas = mixture_model.covariances_.flatten()

    sparsify = ~isclose(ws,0,atol=1e-1) # non-zero components
    ws,mus,sigmas = ws[sparsify],mus[sparsify],sigmas[sparsify]

    sort_index = argsort(mus) # order by means
    ws,mus,sigmas = ws[sort_index],mus[sort_index],sigmas[sort_index]

    # locate extrema
    minima, = argrelextrema(px,less)
    maxima, = argrelextrema(px,greater)

    minima = sort(x[minima])
    maxima = sort(x[maxima])

    # unimodal distributions
    if len(maxima) == 1 :

        # negative populations
        index = argmax(ws)
        if maxima[0]<cutoff :
            threshold = mus[index]+5*sigmas[index]

        # positive populations
        if maxima[0]>cutoff :
            threshold = mus[index]-5*sigmas[index]

        # figure(figsize=(6,6))

        # hist,bins = histogram(channel,bins=200,density=True)
        # bins = (bins[1:]+bins[:-1])/2
        # plot(bins,hist,color='gray',linewidth=1)

        # # negative populations
        # if all(mus<cutoff):

        #     i = 0
        #     for w,m,s in zip(ws,mus,sigmas):

        #         fill_between(x,w*exp(-(x-m)**2/(2*s))/sqrt(2*pi*s),
        #             color='red' if i==0 else 'gray',alpha=0.5)
        #         i += 1


        # # positive populations
        # elif all(mus>cutoff):

        #     i = 0
        #     for w,m,s in zip(ws[::-1],mus[::-1],sigmas[::-1]):

        #         fill_between(x,w*exp(-(x-m)**2/(2*s))/sqrt(2*pi*s),
        #             color='green' if i==0 else 'gray',alpha=0.5)
        #         i += 1

        # elif len(mus) == 2 :

        #     i = 0
        #     for w,m,s in zip(ws,mus,sigmas):
        #         fill_between(x,w*exp(-(x-m)**2/(2*s))/sqrt(2*pi*s),
        #             color='red' if i==0 else 'green',alpha=0.5)
        #         i += 1

        # else :

        #     i = 0
        #     for w,m,s in zip(ws,mus,sigmas):
        #         if i == 0 :
        #             fill_between(x,w*exp(-(x-m)**2/(2*s))/sqrt(2*pi*s),color='red',alpha=0.5)
        #         elif i == 1 :
        #             fill_between(x,w*exp(-(x-m)**2/(2*s))/sqrt(2*pi*s),color='gray',alpha=0.5)
        #         else :
        #             fill_between(x,w*exp(-(x-m)**2/(2*s))/sqrt(2*pi*s),color='green',alpha=0.5)
        #         i += 1

        # axvline(cutoff,color='k')
        # show()


    # bimodal distributions
    if len(maxima) == 2 :
        threshold = minima[0]

    # trimodal distributions
    if len(maxima) == 3 :
        threshold = maxima[1]

    # use stricter threshold 
    if cutoff < threshold :
        return threshold
    else :
        return cutoff