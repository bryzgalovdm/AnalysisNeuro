function [Mean, STD] = GetNanDescStats(dat, dim)
% 
% Wraps both nanmean and nanstd in one function
% 
% INPUT
% 
%     dat        numerical array
%     dim        an integer specifying which dimensions to get mean and std
%                across
% 
%  OUTPUT
%
%     Mean      nanmean of dat across dim's dimension
%     STD       nanstd of dat across dim's dimension
%
% by Dima Bryzgalov from MOBS team, Paris, France
% 29/05/2020
% github.com/bryzgalovdm

if nargin == 1
    dim = 1;
end

Mean = nanmean(dat,dim);
STD = nanstd(dat,[],dim);

end