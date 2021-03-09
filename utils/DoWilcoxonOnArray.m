function [pvalues, pairs_ToCompare] =  DoWilcoxonOnArray(data, pairs_ToCompare)
% 
% Calculates p-values for Wilcoxon ranksum test on selected parts of dataset
% it works with N*M matrices where N is number of observations and M is
% number of samples, or with cell arrays
% 
% INPUT
% 
%    data               data to compare statistically: either N*M matrix where
%                       colums are for comparison, ot cell array where all cells to
%                       compare are of equal size
%    pairs_ToCompare    cell array with pairs to be do the test on. If data
%                       is matrix, columns will be compared, if data is
%                       cell array, individual cells will be compared
% 
% OUTPUT
% 
%    pvalues            p-values for the comprison of each pair
% 
%    pairs_ToCompare    copy of an argument (see above)
% 
% EXAMPLE
% 
%    [pvalues, pairs_ToCompare] =  DoWilcoxonOnArray(data, {[1 2], [1 3], [2 3]});
% 
% Coded by Dima Bryzgalov, MOBS team, Paris,
% 05/03/2021
% github.com/bryzgalovdm
% 

%% Sanaity check on arguments
if ~isa(data, 'cell')
    if ~isa(data, 'numeric')
        error('Incorrect data type (type ''help DoWilcoxonOnArray'' for details).')
    else
        if length(size(a)) ~= 2 && size(a,2) ~= 2
            error('Incorrect data type (type ''help DoWilcoxonOnArray'' for details).')
        else
            celldata = false;
        end
    end
else
    celldata = true;
end

if ~isa(pairs_ToCompare, 'cell')
    error('Incorrect pairs_ToCompare (type ''help DoWilcoxonOnArray'' for details).')
end

%% Do Wilcoxon tests
pvalues = nan(length(pairs_ToCompare), 1);
for ipair = 1:length(pairs_ToCompare)
    if celldata
        pvalues(ipair) = ranksum(data{pairs_ToCompare{ipair}(1)}, data{pairs_ToCompare{ipair}(2)});
    else
        pvalues(ipair) = ranksum(data(:,pairs_ToCompare{ipair}(1)), data(:, pairs_ToCompare{ipair}(2)));
    end
end