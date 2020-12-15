function IsOverlap = FindOverlapArea_PC(PlaceFieldStats, MaskedMap, OverlapFactor)
% 
% Function decided whether given place field overlap with the mask
% 
%  INPUT
% 
%         PlaceFieldStats           Output stats from PlaceField_DB.m
%         MaskedMap                 Mask of the area in question
%         OverlapFactor             Number of pixels overlapping between <MaskedMap> and
%                                   place field that would allow to call place cell
%                                   overlapping with the <MaskedMap> :)
% 
%  OUTPUT
%  
%         IsOverlap                 boolean: true - place field overlaps
%                                   with the area or not
% 
% 
%  EXAMPLE
% 
%           IsOverlap = FindOverlapArea_PC(stats, MaskedMap, 5);
% 
% Coded by Dima Bryzgalov, MOBS team, Paris, France
% 14/12/2020
% github.com/bryzgalovdm

%% Find all overlapping cells
if iscell(PlaceFieldStats.field)
    for k=1:2 % Sometimes there are two place fields
        OverlappedFields = MaskedMap & PlaceFieldStats.field{k};
        numOverlap = nnz(OverlappedFields);
        if numOverlap > OverlapFactor
            IsOverlap = true;
            break
        end
        IsOverlap = false;
    end
else
    OverlappedFields = MaskedMap & PlaceFieldStats.field;
    numOverlap = nnz(OverlappedFields);
    if numOverlap > OverlapFactor
        IsOverlap = true;
    else
        IsOverlap = false;
    end
end

end