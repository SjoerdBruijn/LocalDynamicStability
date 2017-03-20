function [state]=makestatelocal(signal,hc,n_dim)
% build a state space for the calculation of local dynamic stability. a
% standard delay of 10 samples is chosen, and the signal is resampled so
% that on average each stride is 100 samples in length ( see England &
% granata) 
% Input:    signal: the signal you want to make a state space for.
%           hc: vector with moments (in samples) in which heelstrikes happen.
%           ideally, has the same number of elements for each subject/
%           condition, as the number of strides has been shown to influence
%           lds calculations. 
%           n_dim: number of dimensions that the state space should consits of. 
% Output    state: n_dim dimensional state space, where all the states have
%           been normalized such that there length is (length(hc)-1)*100 samples
%           (i.e. 100 samples per stride on average.

%% time normalization
n_strides   = length(hc);
signal_new  = signal(hc(1):hc(end));% take out only the time period we need
t_new       = 0:1:((length(signal_new)-1)); % fix time axis for time normalization
t_interp    = (1:(n_strides-1)*100)/((n_strides-1)*100)*t_new(end); % new time axis, contains (n_steps-1)*100 data points
signal_interp = interp1(t_new,signal_new,t_interp,'spline'); % interpolate
%% create state space
delay = 10;
state = [];
for dim = 1:n_dim
    state = [state signal_interp(((dim-1)*delay)+1:end-(n_dim-dim)*delay)'];
end