function outvector = flatten_cellarray(cellarray, depth)
% 
% Function flattens cell array of vectors or nested cell array of vector
% into one vector
% 
%  INPUT
% 
%         cellarray       array of cell - maximum depth is 2
%         depth           actual depth of the array
%         
%  OUTPUT
%  
%         outvector       vector of all numbers from cell array
% 
%  EXAMPLE
% 
%         cellarray = {{[0 1] [0]}, {[1 2]}, {[2 3 4] [4 5 6]}, {[]}};
%         outvector = flatten_cellarray(cellarray, 2)
% 
% Coded by Dima Bryzgalov, MOBS team, Paris, France
% 03/12/2020
% github.com/bryzgalovdm


%% Flatten
if depth == 1
    if size(cellarray,2) > size(cellarray,1)
        outvector = [cellarray{:}];
    elseif size(cellarray,2) < size(cellarray,1)
        outvector = vertcat(cellarray{:});
    end
elseif depth == 2
    if size(cellarray,2) > size(cellarray,1)
        temp = [cellarray{:}];
    elseif size(cellarray,2) < size(cellarray,1)
        temp = vertcat(cellarray{:});
    end
    
    if size(temp,2) > size(temp,1)
        outvector = [temp{:}];
    elseif size(temp,2) < size(temp,1)
        outvector = vertcat(temp{:});
    end
else
    error('Incorrect value for argument ''depth'': it accepts either 1 or 2');
end

end