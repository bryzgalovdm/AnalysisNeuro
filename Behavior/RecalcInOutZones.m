function [Zone, ZoneLabels, ZoneIndices] = RecalcInOutZones(mask, Xtsd, Ytsd, Ratio_IMAonREAL, varargin)
%
% Calculates zones, one of which go along the borders of the mask and
% another would represent the central part
%
%  INPUT
%
%       mask                mask of the environment
%       Xtsd                tsd with X coordinate
%       Ytsd                tsd with Y coordinate
%       Ratio_IMAonREAL     number that allows to tranlate pixels to cm (comes
%                           from tracking)
%                           Dropbox (boolean - default = false)
%       Ratio               Ration between walls and center (default=0.33)
%
%
%  OUTPUT
%
%       Figure
%
%       See
%
%       UMazeTracking_IRComp_ObjectOrientated_DBparam, PAG_MFB_alongWalls_fig
%
% Coded by Dima Bryzgalov, MOBS team, Paris, France
% 05/02/2021
% github.com/bryzgalovdm

%% Parameters management
Ratio = 0.33;
% Optional parameters
for i=1:2:length(varargin)
    switch(lower(varargin{i}))
        case 'ratio'
            Ratio = varargin{i+1};
            if ~isnumeric(Ratio)
                if Ratio < 0 || Ratio > 1
                    error('Incorrect value for property ''Ratio'' (type ''help RecalcInOutZones'' for details).');
                end
            end
    end
end

%% Creates zones (masks)
stats = regionprops(mask, 'Area');
tempmask=mask;
AimArea=stats.Area*Ratio;
ActArea=stats.Area;
while AimArea<ActArea
    tempmask=imerode(tempmask,strel('disk',1));
    stats = regionprops(tempmask, 'Area');
    ActArea=stats.Area;
end
Zone{1}=uint8(tempmask); % In
Zone{2}=uint8(mask-tempmask);% Out
ZoneLabels={'Inside','Outside'};

%% Calculate ZoneIndices
ZoneIndices = cell(length(Zone), 1);
for izone = 1:length(Zone)
    ZoneIndices{izone}=find(diag(Zone{izone}(floor(Data(Xtsd)*Ratio_IMAonREAL),floor(Data(Ytsd)*Ratio_IMAonREAL))));
end

end