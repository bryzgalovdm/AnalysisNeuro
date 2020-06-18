function [xlimits, ylimits] = DefineGoodFigLimits_2D(X,Y,precise)
% 
% Defines comfortable limits of 2D figure based on limits points of data
% 
% INPUT
% 
%     X             First set of data
%     Y             Second set of data
%     
% OUTPUT
% 
%     xlimits       two elements vector with xlim
%     ylimits       two elements vector with ylim
%     precise(opt)  true if you want limits that correspond to limits of your
%                   data (default - false)
% 
% Coded by Dima Bryzgalov, MOBS team, Paris, France
% 18/06/2020
% github.com/bryzgalovdm
    
%% Argument management and defaults
if nargin < 2
    error('Not enough arguments');
elseif nargin == 2
    precise = false;
elseif nargin == 3
    if ~islogical(precise)
        error('"Precise" should be logical value');
    end
else
    error('Too many arguments');
end

%% The body
% X
minx = min(X);
maxx = max(X);
difx = maxx-minx;
if precise
    xlimits = [minx-difx maxx+difx];
else
    xlimits = [minx-difx*0.05 maxx+difx*0.05];
end

%Y
miny = min(Y);
maxy = max(Y);
dify = maxy-miny;
if precise
    ylimits = [miny-dify*0.05 maxy+dify];
else
    ylimits = [miny-dify*0.05 maxy+dify*0.05];
end

end