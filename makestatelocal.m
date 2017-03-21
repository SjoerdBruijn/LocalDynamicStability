function [state]=makestatelocal(signal,hc,n_dim,delay)
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
%           delay: (optional!) delay in samples. When no input is given,
%           this is set to 1. 
% Output    state: n_dim dimensional state space, where all the states have
%           been normalized such that there length is (length(hc)-1)*100 samples
%           (i.e. 100 samples per stride on average.
%% Copyright
%     COPYRIGHT (c) 2017 Sjoerd Bruijn, VU University Amsterdam
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%% version history; march 2017, SMB updated from older versions by SMB, Jaap van Dieen, Kim van Schooten & others to make suitable
% for general use (i.e. proper commenting, cleanup of code, etc)
%% check for inputs
if nargin <= 3
    delay = 10;
end
n_strides   = length(hc); % how much heelstrikes do we have
n_samples   = (n_strides-1)*100;% amount of samples we will normalize to

%% set up outputs
state = NaN*ones(n_samples-delay*(n_dim-1),n_dim);

%% time normalization
signal_new  = signal(hc(1):hc(end));% take out only the time period we need
t_new       = 1:length(signal_new); % fix time axis for time normalization
t_interp    = (1:n_samples)/n_samples*t_new(end); % new time axis, contains (n_steps-1)*100 data points
signal_interp = interp1(t_new,signal_new,t_interp,'spline'); % interpolate

%% create state space
for i_dim = 1:n_dim
    state(:,i_dim) = signal_interp(((i_dim-1)*delay)+1:end-(n_dim-i_dim)*delay);
end