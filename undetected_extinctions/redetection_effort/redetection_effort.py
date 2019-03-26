import numpy as np
from scipy import interpolate

# calculates the negative log-likelihood
def negLogLikelihoodSpline(cs, ts, spps):
    '''
    cs: list of c_t values defining the spline
    ts: list of t values defining the spline
    spps: dictionary containing species' redetection records
    '''

    eps = 1e-6 # a small number used in place of 0 or 1 probability of redetection

    # create the timeseries and a dictionary for looking up
    t0 = ts[0]; tf = ts[-1]; tV = list(range(t0,tf+1))
    t2i = { t: i for i, t in enumerate( range(t0, tf+2) ) }

    # create the linear spline function using the points defined in ts and cs
    f = interpolate.interp1d( ts, cs, kind='linear')
    cV_fit = f(tV)

    # get the redetection probability estimate for each species
    redetn_probs = { name: min( 1-eps, len(spp['redetns']) / sum( cV_fit[ t2i[max(t0,spp['frst']+1)] : t2i[min(tf+1,spp['last'])] ] ) ) for name, spp in spps.items() }

    # now calculate the likelihood of the data given the c and intrinsic redetection probabilities

    logLtV = list() # store log of probabilities at each timestep
    for t, c in zip(tV, cV_fit):

        # the probability of each timestep of data being how it was if c is true
        pp = list()
        for name, spp in spps.items():

            if t > spp['frst'] and t < spp['last']: # only count species who were available to be redetected

                redetn_prob = redetn_probs[name]

                if t in spp['redetns']:

                    pp.append( redetn_prob*c )

                else:

                    pp.append( 1-redetn_prob*c )

        logLt = sum(np.log(pp)) # rather than using Lt = reduce(operator.mul, pp, 1), to avoid numerical rounding errors

        logLtV.append(logLt)

    negLL = -sum(logLtV)

    return negLL
