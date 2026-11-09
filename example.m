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
tic
[divergence,lds] = lds_calc(state,ws,fs,period, showplot)
toc

