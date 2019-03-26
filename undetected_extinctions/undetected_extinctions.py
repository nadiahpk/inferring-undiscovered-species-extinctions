# 

import numpy as np
from scipy.stats import hypergeom
from itertools import compress
from scipy.special import logit


def frst_last_changed_E(frst_last):
    """
    years_mod, frst_last_mod = first_last_changed_E(frst_last):

    Accepts a list of first and last detection years for species,
    finds the years in which at least one species had its last detection,
    and modifies the years of first detection so that they match the modified list of years.

    frst_last:
        list of tuples of two integers, a list of species first and last detection years
        e.g. [ (1891, 1894), (1895, 2015) ] means the second species was first detected in 1895 and last detected in 2015
    years_mod:
        list of integers, years in which at least one species had its last detection
    frst_last_mod:
        list of tuples of two integers, a list of species first and last detection years where first detection year is modified
        to match a year included in years_mod

    >>> frst_last = [ (0,0), (0,3), (0,6), (2,3), (1,7) ]
    >>> years_mod, frst_last_mod = frst_last_changed_E(frst_last)
    >>> years_mod
    [0, 3, 6, 7]
    >>> frst_last_mod
    [(0, 0), (0, 3), (0, 6), (3, 3), (3, 7)]
    """

    frstObs, lastObs = zip(* frst_last )

    y0 = min(frstObs)
    yf = max(lastObs)
    assert(y0 <= yf) # check we received a sensible input
    years = list(range(y0, yf+1)) # y0, y0+1, ..., yf inclusive

    # years containing at least one last observation

    years_mod = [ l for l in sorted(set(lastObs)) if l < yf ]

    # make sure first and last year is in our set as well

    if y0 not in years_mod:

        years_mod = [y0] + years_mod

    if yf not in years_mod:

        years_mod = years_mod + [yf]

    # move each first observation to the new year where needed

    frstObs_ = list()

    for f in frstObs:

        if f in years_mod:

            frstObs_.append(f)

        else:

            ff = next(y for y in years_mod if f < y)

            frstObs_.append(ff)

    frstObs_ = np.array(frstObs_)

    frst_last_mod = list(zip(frstObs_, lastObs))

    return years_mod, frst_last_mod

def SE_changed_E(S_orig, E_orig):
    '''
    S, E = SE_changed_E(S_orig, E_orig)

    Accepts S, E, vectors (e.g. from a simulation), and modifies the record so that timesteps include
    at least one detected species going extinct. Has a similar purpose to function frst_last_changed_E.

    S_orig:
        numpy array of integers, the total number of detected extant species in each year
    E_orig:
        numpy array of integers, the total number of detected extinct species in each year

    S, E:
        same as above, but with timesteps collapsed so that each timestep includes at least one detected extinction

    >>> S_orig = np.array([3,3,4,4,2,2,2,1])
    >>> E_orig = np.array([0,1,1,1,3,3,3,4])
    >>> S, E = SE_changed_E(S_orig, E_orig)
    >>> S
    array([3, 4, 2, 1])
    >>> E
    array([0, 1, 3, 4])
    '''

    # create the list of indices (years) had at least one detected extinction

    idxs = list(range(len(S_orig))) # list of all indexes
    extn_flags = E_orig[1:] - E_orig[:-1] > 0 # flag idxs in which an extinction occurred
    idxs_mod = list(compress(idxs, extn_flags)) # list of idxs with extinctions

    # put the first and last onto the record if needed

    if idxs[0] not in idxs_mod:
        idxs_mod = [idxs[0]] + idxs_mod
    if idxs[-1] not in idxs_mod:
        idxs_mod = idxs_mod + [idxs[-1]]

    # return only those years

    S = S_orig[idxs_mod]; E = E_orig[idxs_mod];

    return S, E

def get_SE(frst_last, years=None):
    """
    years, S, E = get_SE(frst_last, years=None)

    Accepts a list of first and last observations for species,
    and optionally a list of years,
    and calculate the number of extant (S) and extinct (E) species in timeseries.

    frst_last:
        list of tuples of two integers, a list of species first and last detection years
        e.g. [ (1891, 1894), (1895, 2015) ] means the second species was first detected in 1895 and last detected in 2015
    years:
        list of years, corresponding to S and E
    S:
        numpy array of integers, the total number of detected extant species in each year
    E:
        numpy array of integers, the total number of detected extinct species in each year
    """

    # separate the first and last observations and turn into numpy arrays

    frstObs, lastObs = zip(* frst_last )
    frstObs = np.array(frstObs)
    lastObs = np.array(lastObs)

    # sort out the range of years for which we create the timeseries

    if years is None: # use every year in the range

        y0 = min(frstObs); yf = max(lastObs);
        assert(y0 <= yf) # check we received a sensible input
        years = list(range(y0, yf+1)) # y0, y0+1, ..., yf inclusive

    # calculate S and E for each year in timeseries

    S = np.array( [ sum( (frstObs <= t) & (lastObs >= t) ) for t in years ] )
    E = np.array( [ sum( lastObs < t ) for t in years ] )

    return years, S, E

def find_U0_bnd(alpha, S0, S1, U1, d0, impossibleFlag=False, omega=None, biasedurn=None):
    '''
    U0_bnd = inverse_midp(alpha, min_poss_U0, S0, S1, U1, d0, omega=None, biasedurn=None)

    Find the bound for the number of undetected extant species at the previous timestep (U0)
    given the data at a given timestep (S1, U1, d0 etc.) and confidence level using the
    mid-P method.

    By default the model used is the central hypergeometric distribution. When omega is specified,
    the Fisher non-central hypergeometric distribution is used instead.

    alpha:
        float, confidence level (e.g. obtained by randomly sampling ~ U(0,1))
    S0:
        integer, the number of detected extant species at the previous timestep
    S1:
        integer, the number of detected extant species at the current timestep
    U1:
        integer, the number of undetected extant species at the current timestep
    d0:
        integer, the number of species detected during the previous timestep
    impossibleFlag:
        logical, whether or not our data is already in the "impossible" region
    omega:
        float, the odds ratio of survival in undetected / detected species
    biasedurn:
        rpy2.robjects.packages.Package as a <module 'BiasedUrn'>, used to access
        the R package BiasedUrn via the rpy2 package so that the Fisher variant
        can be modelled. It needs to be loaded and passed from the function that
        calls this function, allong with the objects passer, 
        i.e. rpy2.robjects.numpy2ri.activate(); biasedurn = rpy2.robjects.packages.importr('BiasedUrn')
    '''

    min_poss_U0 = U1 + d0 # the minimum possible value of U0 in reality


    # define the mid-P function that we will be using
    # ---

    if omega: # doing the Fisher variant

        midP_fnc = lambda U0, S0, S1, U1, d0: 0.5 * ( biasedurn.pFNCHypergeo( U1+d0, U0, S0, S1+U1, omega )[0] + biasedurn.pFNCHypergeo( U1+d0-1, U0, S0, S1+U1, omega )[0] )

    else: # central hypergeom, assumes equal probability of extinction

        midP_fnc = lambda U0, S0, S1, U1, d0: 0.5 * ( hypergeom.cdf( U1+d0, S0+U0, U0, S1+U1 ) + hypergeom.cdf( U1+d0-1, S0+U0, U0, S1+U1 ) )



    # obtain a sample value of U0 at our confidence level alpha
    # ---

    # first, check if we're in the situation where we wouldn't accept the minimum possible value of U0

    if midP_fnc(min_poss_U0, S0, S1, U1, d0) < alpha:

        if impossibleFlag:

            # don't take two steps in to the impossible region
            U0_bnd = min_poss_U0

        else:

            # U0 is set to the 'impossible' value, one less than minimum value, end search
            U0_bnd = min_poss_U0-1 
            impossibleFlag = True

    else:

        # We are trying to find the greatest value of U0, called U0_bnd, such that
        # midP_fnc(U0_bnd, ...) > alpha. 
        # In other words, we need to find the value such that the mid-P function switches:
        # midP_fnc(U0_bnd, ...) > alpha  and  midP_fnc(U0_bnd+1, ...) < alpha.
        # The search has two phases:
        # 1. find an interval [U0_lo, U0_hi] within which U0_bnd lies
        # 2. binary search between our upper and lower search bounds, U0_lo and U0_hi, to find U0_bnd


        # 1. find an interval [U0_lo, U0_hi] within which U0_bnd lies
        # ----

        # To quickly find the interval, we start at the minimum possible value for U0_bnd,
        # and work our way upwards until we find a U0_hi such that midP_fnc(U0_hi, ...) < alpha.
        # The step-size from U0_lo to the next U0 value is doubled each time, to account
        # for the fact that there is no upper limit to the U0 value

        # initialise for the loop
        U0_lo = min_poss_U0
        step_size = 1 # this will double at each step
        U0_hi_found = False

        while not U0_hi_found:

            # prospective upper
            U0_hi = U0_lo + step_size

            if midP_fnc(U0_hi, S0, S1, U1, d0) < alpha:

                U0_hi_found = True  # we've found an upper bound to search within

            else:

                U0_lo = U0_hi   # we didn't find an upper so this is new lower bound for search
                step_size *= 2  # double the step size for next time


        # 2. binary search between our upper and lower search bounds, U0_lo and U0_hi
        # ----

        # Once an interval [U0_lo, U0_hi] has been identified, we use a binary search
        # to find the value of U0_lo such that 
        # midP_fnc(U0_lo, ...) > alpha  and  midP_fnc(U0_lo+1, ...) < alpha.
        # This value of U0_lo is U0_bnd.
        # The binary search works by iteratively taking the mid-point U0_mid between U0_lo and U0_hi,
        # checking if it is a new low or high bound on the interval, and checking the condition above
        # to terminate the search.

        # have we already found bound?
        U0_bnd_found = U0_hi - U0_lo == 1 # if they're one apart, then we've found the switch point

        while not U0_bnd_found:

            U0_mid = int( (U0_lo+U0_hi) / 2 ) # floored midpoint
            alpha_mid = midP_fnc(U0_mid, S0, S1, U1, d0) # the alpha value at this mid-point

            if alpha_mid == alpha: # it's the actual bound (unlikely to happen)

                U0_bnd_found = True
                U0_lo = U0_mid # the bound is stored in U0_lo

            else:

                if alpha_mid < alpha: # the mid-point is a new upper bound

                    U0_hi = U0_mid

                else: # elif alpha_mid > alpha: # the mid-point is a new lower bound

                    U0_lo = U0_mid

                # if the interval bounds are one apart, then we've found U0_bnd
                U0_bnd_found = U0_hi - U0_lo == 1 

        # When we exit the while loop above then the bound has been found.
        # It is stored in U0_lo (i.e. the greatest value of U0 for which midP_fnc >= alpha)
        U0_bnd = U0_lo

        if U0_bnd > min_poss_U0:
            # we've moved out of the impossible region
            impossibleFlag = False

    if U0_bnd < 0:
        U0_bnd = 0 # don't allow negative numbers of undetected species



    return U0_bnd, impossibleFlag
    
# old version of find_U0_bnd, use find_U0_bnd instead
def inverse_midp(alpha, min_poss_U0, S0, S1, U1, d0, omega=None, biasedurn=None):
    '''
    U0_bnd = inverse_midp(alpha, min_poss_U0, S0, S1, U1, d0, omega=None, biasedurn=None)

    Find the bound for the number of undetected extant species at the previous timestep (U0)
    given the data at a given timestep (S1, U1, d0 etc.) and confidence level using the
    mid-P method.

    By default the model used is the central hypergeometric distribution. When omega is specified,
    the Fisher non-central hypergeometric distribution is used instead.

    alpha:
        float, confidence level (e.g. obtained by randomly sampling ~ U(0,1))
    min_poss_U0:
        integer, the minimum possible value of undetected extant species at the previous timestep
    S0:
        integer, the number of detected extant species at the previous timestep
    S1:
        integer, the number of detected extant species at the current timestep
    U1:
        integer, the number of undetected extant species at the current timestep
    d0:
        integer, the number of species detected during the previous timestep
    omega:
        float, the odds ratio of survival in undetected / detected species
    biasedurn:
        rpy2.robjects.packages.Package as a <module 'BiasedUrn'>, used to access
        the R package BiasedUrn via the rpy2 package so that the Fisher variant
        can be modelled. It needs to be loaded and passed from the function that
        calls this function, allong with the objects passer, 
        i.e. rpy2.robjects.numpy2ri.activate(); biasedurn = rpy2.robjects.packages.importr('BiasedUrn')
    '''


    # define the mid-P function that we will be using
    # ---

    if omega: # doing the Fisher variant

        midP_fnc = lambda U0, S0, S1, U1, d0: 0.5 * ( biasedurn.pFNCHypergeo( U1+d0, U0, S0, S1+U1, omega )[0] + biasedurn.pFNCHypergeo( U1+d0-1, U0, S0, S1+U1, omega )[0] )

    else: # central hypergeom, assumes equal probability of extinction

        midP_fnc = lambda U0, S0, S1, U1, d0: 0.5 * ( hypergeom.cdf( U1+d0, S0+U0, U0, S1+U1 ) + hypergeom.cdf( U1+d0-1, S0+U0, U0, S1+U1 ) )


    # obtain a sample value of U0 at our confidence level alpha
    # ---

    # first, check if we're in the situation where we wouldn't accept the minimum possible value of U0

    if midP_fnc(min_poss_U0, S0, S1, U1, d0) < alpha:

        U0_bnd = min_poss_U0-1 # U0 is set to the 'impossible' value, one less than minimum value, end search

    else: # start at the minimum possible value

        # the function we'll use here on out, logit because slope changes less and it's unbounded
        logit_alpha = logit(alpha)
        f = lambda U0: logit( midP_fnc(U0, S0, S1, U1, d0) ) - logit_alpha

        # initialise for my loop
        U0_next = min_poss_U0
        f_next = f(U0_next)

        # For small values of omega, the mid-P function has rounding issues and may return f_next = infinity. 
        # Can't just step forward until that's not true (i.e. while f_next == np.inf: ), as sometimes 
        # returns back to inf after returning finite value. So set a minimum value on f
        while f_next > 20: # ~ logit( 1 - 1e-9 )
            U0_next += 1; f_next = f(U0_next)

        notfound_U0_bnd = True

        # loop using the Newton's method to find sample value
        while notfound_U0_bnd:

            # update position

            U0_prev = U0_next; f_prev = f_next

            # do a Newton's step

            df = f(U0_prev+1) - f_prev              # get slope at new position
            U0_next = int( U0_prev - f_prev / df )  # next position
            f_next = f(U0_next)                     # value at next position

            # check if I've found it

            f_next_befor = f(U0_next-1); f_next_after = f(U0_next+1);

            if np.sign(f_next_befor) != np.sign(f_next):

                U0_bnd = U0_next-1 # always the smaller because this function decreases w increasing U0
                notfound_U0_bnd = False
        
            if np.sign(f_next) != np.sign(f_next_after):

                U0_bnd = U0_next
                notfound_U0_bnd = False

    if U0_bnd < 0:
        U0_bnd = 0 # don't allow negative numbers of undetected species

    return U0_bnd
