function [opairs, dpairs, udpairs] = FindOverlappingPlaceFields(S, Xtsd, Ytsd, epoch, overlapFactor, varargin)
%%% FindOverlappingPlaceFields
%
% Splits all the possible pairs of place cells into opairs (place cells with
% overlapping place fields) and dpairs (place cells with non-overlapping place fields)
%
% INPUT
%
%     S              tsd array of spike train (with place cells - for good performance)
%     Xtsd           x-coordinates for place field definition (in tsd)
%     Ytsd           y-coordinates for place field definition (in tsd)
%     epoch          epoch in which you would like to assess overlapping (intervalSet)
%     overlapFactor  how many pixels in place fields should overlap to become <opairs>
%
%     optional arguments:
%
%     smoothing      spatial smooting factor (default - 2)
%     sizemap        size of the map (it will be size * size) in pixels (default - 50)
%
%  OUTPUT
%
%     opairs         cell consisting of 2-numbers vectors with indices of
%                    overlapping pairs
%     dpairs         cell consisting of 2-numbers vectors with indices of
%                    non-overlapping (distant) pairs
%
%  NOTE
%
%    Please use spike trains of place cells only to ensure good performance
%
%  SEE ALSO
%
%   PlaceField_DB   FindAllPlaceCells
%
% By Dima Bryzgalov, MOBS team, Paris,
% 30/04/2020
% github.com/bryzgalovdm

%% Optional arguments
for i=1:2:length(varargin)
    
    switch(lower(varargin{i}))
        
        case 'smoothing'
            smoothing = varargin{i+1};
            if ~isa(smoothing,'numeric')
                error('Incorrect value for property ''Smoothing'' (type ''help PlaceField'' for details).');
            end
        case 'sizemap'
            sizeMap = varargin{i+1};
            if ~isa(sizeMap,'numeric')
                error('Incorrect value for property ''Size'' (type ''help PlaceField_DB'' for details).');
            end
            
    end
end

% Default values of optional arguments
try
    smoothing;
catch
    smoothing=2;
end

try
    sizeMap;
catch
    sizeMap=50;
end

%% PreAllocation
stats = cell(1,length(S));
num_possile_pairs = ceil(factorial(length(S))/(factorial(2)*factorial(length(S)-2)));
opairs = cell(1, num_possile_pairs);
dpairs = cell(1, num_possile_pairs);
udpairs = cell(1, num_possile_pairs);

%% Find the place fields
for i = 1:length(S)
    try
        [~, ~, stats{i}]=PlaceField_DB(S{i}, Xtsd, Ytsd,...
            'Epoch', epoch, 'SizeMap', sizeMap, 'Smoothing', smoothing, 'PlotResults', 0, 'PlotPoisson', 0);
    catch
        stats{i}=[];
    end
end

%% FindOverlappingPlaceCells
cnto = 1;
cntd = 1;
cntu = 1;
for j=1:length(S)
    for k = j+1:length(S)
        if isempty(stats{j}) || isempty(stats{j}.field) || isempty(stats{k}) || isempty(stats{k}.field)
            udpairs{cntu} = [j k];
            cntu = cntu + 1;
        else
            if iscell(stats{j}.field) && ~iscell(stats{k}.field) % If cell1 has two fields, and cell2 has one
                OverlappedFields{1} = stats{j}.field{1} & stats{k}.field;
                OverlappedFields{2} = stats{j}.field{2} & stats{k}.field;
                numOverlap{1} = nnz(OverlappedFields{1});
                numOverlap{2} = nnz(OverlappedFields{2});
                % If either one ot the other place field overlaps, consider overlap
                if numOverlap{1} > overlapFactor || numOverlap{2} > overlapFactor
                    opairs{cnto} = [j k];
                    cnto=cnto+1;
                else
                    dpairs{cntd} = [j k];
                    cntd=cntd+1;
                end
            elseif ~iscell(stats{j}.field) && iscell(stats{k}.field) % If cell1 has one field, and cell2 has two
                OverlappedFields{1} = stats{j}.field & stats{k}.field{1};
                OverlappedFields{2} = stats{j}.field & stats{k}.field{2};
                numOverlap{1} = nnz(OverlappedFields{1});
                numOverlap{2} = nnz(OverlappedFields{2});
                if numOverlap{1} > overlapFactor || numOverlap{2} > overlapFactor
                    opairs{cnto} = [j k];
                    cnto=cnto+1;
                else
                    dpairs{cntd} = [j k];
                    cntd=cntd+1;
                end
            elseif iscell(stats{j}.field) && iscell(stats{k}.field) % If both cell1 and cell2 have two fields
                OverlappedFields{1} = stats{j}.field{1} & stats{k}.field{1};
                OverlappedFields{2} = stats{j}.field{1} & stats{k}.field{2};
                OverlappedFields{3} = stats{j}.field{2} & stats{k}.field{1};
                OverlappedFields{4} = stats{j}.field{2} & stats{k}.field{2};
                numOverlap{1} = nnz(OverlappedFields{1});
                numOverlap{2} = nnz(OverlappedFields{2});
                numOverlap{3} = nnz(OverlappedFields{3});
                numOverlap{4} = nnz(OverlappedFields{4});
                if numOverlap{1} > overlapFactor || numOverlap{2} > overlapFactor ||...
                        numOverlap{3} > overlapFactor || numOverlap{4} > overlapFactor
                    opairs{cnto} = [j k];
                    cnto=cnto+1;
                else
                    dpairs{cntd} = [j k];
                    cntd=cntd+1;
                end
            else
                OverlappedFields = stats{j}.field & stats{k}.field; % If both cell1 and cell2 have one field
                numOverlap = nnz(OverlappedFields);
                if numOverlap > overlapFactor
                    opairs{cnto} = [j k];
                    cnto=cnto+1;
                else
                    dpairs{cntd} = [j k];
                    cntd=cntd+1;
                end
            end
            clear OverlappedFields numOverlap
        end
    end
end

%% Clean the resulting variables
opairs = opairs(~cellfun('isempty',opairs));
dpairs = dpairs(~cellfun('isempty',dpairs));
udpairs = udpairs(~cellfun('isempty',udpairs));

end