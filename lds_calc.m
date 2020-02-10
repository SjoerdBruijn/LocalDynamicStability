function [divergence,lds]=lds_calc(state,ws,fs,period, plotje)
% calculate local dynamic stability ( max lyapunov exponent), according to Rosenstein (1993) algorithm
% Input:    state: appropriate state space
%           ws: window size over which divergence should be calculated(in seconds/cycles)
%           fs: sample frequency
%           period: dominant period in the signal(in seconds/cycles),
%           plotje: 1= show a graph, other= don't show graph
% Output:   divergence: the divergence curve
%           lds: the 2 estimates of the local divergence exponents ( long term and short term)
%           note that if ws is not larger then 4*period, long term is not calculated,
%           as this is normally calculated 4*period to ws.
% Note:     The slope calculated is from 0-0.5 strides, and from 4-10
%           strides. when trying to calculate a 'real' Lyapunov exponent from a
%           system, you will have to calculate the slope at an appropriate point in
%           the divergence curve yourself. 
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

%% setting up some variables
ws          = round(ws*fs); % as ws is input in cycles, and better to have it in samples
[m,n]       = size(state);
state       = [state;NaN*ones(ws,n)]; % we extend the state space with NaN, so that we don't run into problems later
divergence  = NaN*ones(m,ws); % set up the output divergence matrix
difference  = NaN*ones(m+ws,n); % set up some difference matrix that we will use

%% find and track nearest neighbours
for i_t = 1:m % loop over time samples
    for i_d = 1:n % loop over dimensions
        difference(:,i_d) = (state(:,i_d)-state(i_t,i_d)).^2;
    end
    start_index         = round(max([1,i_t-round(0.5*period*fs)])); %find point half a period befor current point
    stop_index          = round(min([m,i_t+round(0.5*period*fs)])); %find point half a period past current point
    difference(start_index:stop_index,:) = NaN;% discard data within one period from sample i_t putting it to nan
    [~,index]           = min(sum(difference,2));% find nearest neighbour
    divergence(i_t,:)   = sqrt(sum((state(i_t:i_t+ws-1,:)-state(index:index+ws-1,:)).^2,2)');% track divergence, and store
end
divergence=nanmean (log(divergence)); % calculate average for output

%% calculate least squares fit
L1 = 0.5*period*fs;
Ps = polyfit(1/fs:1/fs:L1/fs,divergence(:,1:L1),1);
if ws>4*period*fs
    L2 = 4*period*fs;
    Pl = polyfit(L2/fs:1/fs:ws/fs,divergence(:,L2:ws),1);
else
    Pl(:,1)=NaN;
end
lds=[Ps(:,1) Pl(:,1)];

%% plot if indicated
if nargout==0 ||plotje==1
    figure;
    Ys = polyval(Ps,1/fs:1/fs:L1/fs);
    Yl = polyval(Pl,L2/fs:1/fs:ws/fs);
    plot((1:ws)/fs,divergence); hold on
    plot(1/fs:1/fs:L1/fs,Ys,'m');hold on
    plot(L2/fs:1/fs:ws/fs,Yl,'r'); hold on
    legend({'Divergence Curve  ';['Lambda S: ',num2str(Ps(1))];['Lambda L: ',num2str(Pl(1)) ]});
    title('Divergence curve');
    xlabel('Time (sec)');
    ylabel('Ln(divergence)');
    set(gca,'Box','off')
end
