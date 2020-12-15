function ripprob = RipProbDB(spiketrain, ripIS)
% 
% Function calculates probablility of a given neuron to fire during
% ripples, i.e. number of ripples where a neuron fired at least once
% divided by overall number of ripples
% 
%  INPUT
% 
%         spiketrain      spikes in tsd format
%         ripIS           intervalSet with start and end of each ripple
%         
%  OUTPUT
%  
%         ripprob         ripples probability
% 
%  EXAMPLE
% 
%           ripprob = RipProbDB(spiketrain, ripIS)
% 
% Coded by Dima Bryzgalov, MOBS team, Paris, France
% 03/12/2020
% github.com/bryzgalovdm

%%
% Convert back to original matlab arrays
sp_times = Range(spiketrain);
rip_is = [Start(ripIS) End(ripIS)];

% Find how many unique ripples contain at least one spike
[~,interval] = InIntervals(sp_times,rip_is);
ripprob = length(unique(nonzeros(interval)))/length(Start(ripIS));

end