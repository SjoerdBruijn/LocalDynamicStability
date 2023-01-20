# LocalDynamicStability
if you use this code, feel free to cite it's DOI; [![DOI](https://zenodo.org/badge/85585950.svg)](https://zenodo.org/badge/latestdoi/85585950)

Calculate local dynamic stability using Rosenstein's (1993) algorithm, the code used for many of my papers. For why it may be a good gait stability measure, see Bruijn et al., (2013).

In short, the code can be used to generate a state space using embedding delay using [makestatelocal.m](https://github.com/SjoerdBruijn/LocalDynamicStability/blob/master/makestatelocal.m) (while normalizing time, as is common in gait research, see Bruijn et al., (2013)), after which the average logarithmic rate of divergence (i.e. the actual local dynamic stability metric) can be calculated either using [lds_calc.m](https://github.com/SjoerdBruijn/LocalDynamicStability/blob/master/lds_calc.m), which uses the traditional Rosenstein (1993) algorithm, or using [lds_calc_mehdizadeh.m](https://github.com/SjoerdBruijn/LocalDynamicStability/blob/master/lds_calc_mehdizadeh.m), which allows for using multiple nearest neighbours, and may be preferable according to Mehdizadeh (2019).

Example code:
```
clear all;
load 'testdata'

ws      = 10;
fs      = 100;
showplot= 1;
period  = 1;
n_dim   = 5;
delay   = 10;

% we want to do this on MLCoM velocity, not position
CoM_ML_vel  = gradient(CoM_ML,1/fs_opto);

% create a 5 dimensional time normalised state space, with delay of 10
% samples
state       = makestatelocal(CoM_ML_vel,events.lhs,n_dim,delay);
% calculate local divergence exponent
[divergence,lds] = lds_calc(state,ws,fs,period, showplot)
```
References:
  - Bruijn, S.M., O.G. Meijer, P.J. Beek, et al., Journal of the Royal Society, Interface (2013) Assessing the stability of human locomotion: a review of current measures. 10(83): p. 20120999.
  - Rosenstein, M.T., J.J. Collins, and C.J. Deluca, Physica D, (1993) A Practical Method for Calculating Largest Lyapunov Exponents from Small Data Sets. 65(1-2): p. 117-134.
  - Mehdizadeh, S., J Biomech, (2019) A robust method to estimate the largest Lyapunov exponent of noisy signals: A revision to the Rosenstein's algorithm. 85: p. 84-91.
