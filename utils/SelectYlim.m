function ylim_out = SelectYlim(fhandle)
% 
% Outputs ylimiits that would fit all axis on the figure
% Selects minimal of all values towards negative Inf and maximum towards
% positive Inf
% 
% INPUT
% 
%    fhandle    handle for the figure
% 
% OUTPUT
% 
%    ylim_out   ylimits to be applied to all axes of fhandle
% 
% EXAMPLE
% 
%    ylim_out = SelectYlim(fhandle)
% 
% Coded by Dima Bryzgalov, MOBS team, Paris,
% 05/03/2021
% github.com/bryzgalovdm
% 

%%
% Get ylimits of each axis
allAxesInFigure = findall(fhandle,'type','axes');
ylims = nan(length(allAxesInFigure), 2);
for iaxis = 1:length(allAxesInFigure)
    ylims(iaxis,:) = ylim(allAxesInFigure(iaxis));
end

% Select the extreme values
ylim_out(1) = min(ylims(:,1));
ylim_out(2) = max(ylims(:,2));

end