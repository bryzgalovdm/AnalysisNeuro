function [map,mapNS,stats,px,py,FR]=PlaceField_DB(tsa,XS,YS,varargin)

% INPUT
% 
%     tsa            spike times in tsd
%     XS             X position in tsd
%     YS             Y position is tsd
%     
%     optional arguments:
%     
%     smoothing      spatial smooting factor (default - 3)
%     freqVideo      sampling rate of the video (default - 15)
%     threshold      to define place field which is area around peak where FR >= threshold*peak (default - 0.7)
%     size           size of the map (it will be size * size) in pixels
%     LargeMatrix    if 1 edges will be added to maps (default - 1 or True)
%     PlotResults    plot the main figure (occupancy, spike counts, rate maps + movements) (default - 1)
%     PlotPoisson    plot the figure that compares the actual rate map with rate maps of the poisson-distributed set of the same data (default - 1)


%  OUTPUT
%
%     map.rate       average firing rate map (in Hz)
%     map.time       occupancy map (in s)
%     map.count      firing count map
%     mapNS          same three fields non-smoothed
%     stats.x        abscissa of the peak (in bins)
%     stats.y        ordinate of the peak (in bins)
%     stats.peak     in-field peak firing rate (in Hz)
%     stats.mean     in-field mean firing rate (in Hz)
%     stats.size     firing field size (in bins)
%     stats.field    firing field (1 = bin in field, 0 = bin not in field)
%     stats.fieldX   firing field x boundaries (in bins)
%     stats.fieldY   firing field y boundaries (in bins)
%
%  NOTE
%
%    Position (x,y) is used to compute column bin(x) and row bin(y) in all
%    output matrices (e.g. 'map.rate'). Increasing values y are upwards.
%
%
%
% Compute the spatial specificity of the firing map, based on the formula proposed
% by Skaggs et al. (1993).
%
% T = sum(sum(map.time));
% if T == 0,
%   stats.specificity = 0;
% else
%   occupancy = map.time/(T+eps);
%   meanFiringRate = sum(sum(map.count))/(sum(sum(map.time)+eps));
%   if meanFiringRate == 0,
%     stats.specificity = 0;
%   else
%     logArg = map.count/meanFiringRate;
%     logArg(logArg <= 1) = 1;
%     stats.specificity = sum(sum(map.count.*log2(logArg).*occupancy))/meanFiringRate;
%   end
% end
%
% 
% Corrected by Dima Bryzgalov on the basis of common code from MOBS team, Paris, France
% 02/04/2020

%% Parameters handling

for i=1:2:length(varargin)
    
    switch(lower(varargin{i}))
        
        case 'smoothing'
            smo = varargin{i+1};
            if ~isa(smo,'numeric')
                error('Incorrect value for property ''Smoothing'' (type ''help PlaceField'' for details).');
            end
            
        case 'video'
            freqVideo = varargin{i+1};
            if ~isa(freqVideo,'numeric')
                error('Incorrect value for property ''Video'' (type ''help PlaceField_DB'' for details).');
            end
            
        case 'threshold'
            threshold = varargin{i+1};
            if ~isa(threshold,'numeric')
                error('Incorrect value for property ''Threshold'' (type ''help PlaceField_DB'' for details).');
            end
            
        case 'size'
            sizeMap = varargin{i+1};
            if ~isa(sizeMap,'numeric')
                error('Incorrect value for property ''Size'' (type ''help PlaceField_DB'' for details).');
            end
            
        case 'largematrix'
            LMatrix = varargin{i+1};
            if ~isa(sizeMap,'logical') || ~isa(sizeMap,'numerical') 
                error('Incorrect value for property ''LargeMatrix'' (type ''help PlaceField_DB'' for details).');
            end
            
        case 'plotresults'
            plotresults = varargin{i+1};
            if ~(plotresults == 1 || plotresults == 0)
                error('Incorrect value for property ''PlotResults'' (type ''help PlaceField_DB'' for details).');
            end
            
        case 'plotpoisson'
            plotpoisson = varargin{i+1};
            if ~(plotpoisson == 1 || plotpoisson == 0)
                error('Incorrect value for property ''PlotPoisson'' (type ''help PlaceField_DB'' for details).');
            end
            
    end
end

% Default values of optional arguments

try
    smo;
catch
    smo=3;
end

try
    freqVideo;
catch
    freqVideo=15;
end

try
    threshold;
catch
    threshold = 0.7;
end

try
    sizeMap;
catch
    sizeMap=50;   
end

try
    LMatrix;
catch
    LMatrix=true;   
end

try
    plotresults;
catch
    plotresults=1;
end

try
    plotpoisson;
catch
    plotpoisson=1;
end


%% Create 2D histogram of occupancy

px =Data(Restrict(XS,tsa,'align','closest'));
py =Data(Restrict(YS,tsa,'align','closest'));

[occH, x1, x2] = hist2d(Data(XS), Data(YS), sizeMap, sizeMap);

% Restrict movement arrays by spike times
pX = Restrict(XS,tsa,'align','closest');
pY = Restrict(YS, tsa,'align','closest');

%% Create 2D histogram of spike counts
pfH = hist2d(Data(pX), Data(pY), x1, x2);

%% Create rate maps
pf = freqVideo * pfH./occH;

% Add edges to the map
if LMatrix
    largerMatrix = zeros(sizeMap+floor(sizeMap/4),sizeMap+floor(sizeMap/4));
    largerMatrix(1+floor(sizeMap/8):sizeMap+floor(sizeMap/8),1+floor(sizeMap/8):sizeMap+floor(sizeMap/8)) = pf;
    pf = largerMatrix;
end

% Clean rate map from nans
warning on
pf(isnan(pf))=0;
sg = sort(pf(~isnan(pf(:))));
th = sg(end-5);
pf(pf>th) = 0;

mapNS.rate = pf';

% Smooth the map
pf=SmoothDec(pf,[smo,smo]);

pf=pf';
map.rate=pf;

%% Add edges to occupancy map
if LMatrix
    largerMatrix = zeros(sizeMap+floor(sizeMap/4),sizeMap+floor(sizeMap/4));
    largerMatrix(1+floor(sizeMap/8):sizeMap+floor(sizeMap/8),1+floor(sizeMap/8):sizeMap+floor(sizeMap/8)) = occH';
    largerMatrix=SmoothDec(largerMatrix,[smo,smo]);
    mapNS.time = largerMatrix/freqVideo;
    map.time=largerMatrix/freqVideo;
else
    mapNS.time = occH';
    occH=SmoothDec(occH'/freqVideo,[smo,smo]);
    map.time=occH;
end

%% Add edges to spike count map
if LMatrix
    largerMatrix = zeros(sizeMap+floor(sizeMap/4),sizeMap+floor(sizeMap/4));
    largerMatrix(1+floor(sizeMap/8):sizeMap+floor(sizeMap/8),1+floor(sizeMap/8):sizeMap+floor(sizeMap/8)) = pfH';
    map.count = largerMatrix;
    largerMatrix=SmoothDec(largerMatrix,[smo,smo]);
    map.count=largerMatrix;
else
    map.count = pfH';
    pfH=SmoothDec(pfH',[smo,smo]);
    map.count=pfH;
end

%% Update size of the maps
if LMatrix
    nBins =sizeMap+floor(sizeMap/4);
else
    nBins =sizeMap;
end


nBinsX = nBins(1);
if length(nBins) == 1
    nBinsY = nBinsX;
    nBins(2) = nBins;
else
    nBinsY = nBins(2);
end

%% Create statistics structure

stats.x = [];
stats.y = [];
stats.field = [];
stats.size = [];
stats.peak = [];
stats.mean = [];
stats.fieldX = [];
stats.fieldY = [];
stats.spatialInfo=[];
stats.sparsity=[];
stats.specificity = [];

% Compute the spatial specificity of the firing map, based on the formula proposed
% by Skaggs et al. (1993).

T = sum(sum(map.time));
if T == 0
    stats.specificity = 0;
else
    occupancy = map.time/(T+eps);
    meanFiringRate = sum(sum(map.count))/(sum(sum(map.time)+eps));
    if meanFiringRate == 0
        stats.specificity = 0;
    else
        logArg = map.count/meanFiringRate;
        logArg(logArg <= 1) = 1;
        stats.specificity = sum(sum(map.count.*log2(logArg).*occupancy))/meanFiringRate;
    end
end

% Determine the firing field (i.e. the connex area where the firing rates are > threshold*peak).
% There are two ways to do this:
% 1) If the location of the center was provided (see 'center' property), the firing field is
%    found around the center.
% 2) Otherwise, the firing field is assumed to be around the bin with maximal firing rate
%    (when there are two fields, the one with the highest firing rate is selected unless the other
%    one is >20% bigger).

if max(max(map.rate)) == 0
    stats.field = false(size(map.rate));
    return;
end

% First clean the map from very small regions of elevated firing rates
map.rate = CleanMap(nBinsX,nBinsY,threshold,map.rate,20);

% In the absence of any information regarding the location of the center,
% find two candidate firing fields
rate = map.rate;
for i = 1:2
    % Find firing field and compute all parameters
    peak(i) = max(max(rate));
    peakLocation = FindPeakLocation(rate);
    x(i) = peakLocation(1);
    y(i) = peakLocation(2);
    field{i} = FindFiringField(x(i),y(i),nBinsX,nBinsY,threshold*peak,rate);
    fieldSize(i) = sum(sum(field{i}));
    % Remove this firing field from the rate map for next iteration
    rate(logical(field{i})) = 0;
end
% Choose between the two candidate fields
if fieldSize(2) == 0
    winner = 1;
else
    sizeRelativeDifference = fieldSize(2)/fieldSize(1);
    if sizeRelativeDifference > 3 % Used to be 0.2
        winner = 2; % Choose firing field #2 (see below)
    elseif sizeRelativeDifference < 0.33 % used to be else
        winner = 1; % Choose firing field #1 (see below)
    else
        winner = [1 2];
    end
end


% Set stats
if length(winner) == 1
    stats.x = x(winner);
    stats.y = y(winner);
    stats.field = logical(field{winner});
    stats.size = fieldSize(winner);
    stats.peak = peak(winner);
    stats.mean = mean(mean(map.rate(stats.field)));
    i = find(max(stats.field,[],1));
    if ~isempty(i)
        stats.fieldX = [i(1) i(end)];
    end
    i = find(max(stats.field,[],2));
    if ~isempty(i)
        stats.fieldY = [i(1) i(end)];
    end
else
    for i=1:2
        stats.x(i) = x(winner(i));
        stats.y(i) = y(winner(i));
        stats.field{i}=field{i};
        stats.size(i) = fieldSize(winner(i));
        stats.peak(i) = peak(winner(i));
        stats.mean(i) = mean(mean(map.rate(logical(stats.field{i}))));
        j = find(max(stats.field{i},[],1));
        if ~isempty(j)
            stats.fieldX{i} = [j(1) j(end)];
        end
        j = find(max(stats.field{i},[],2));
        if ~isempty(j)
            stats.fieldY{i} = [j(1) j(end)];
        end
    end
end

% Some final parameters
stat=SpatialInfo(map);
stats.spatialInfo=stat.info;
stats.sparsity=stat.sparsity;

if iscell(stats.field)
    for j=1:length(stats.field)
        C{j}=GravityCenter(stats.field{j}.*map.rate);
        GC{j}(j,1)=C{j}(1);
        GC{j}(j,2)=C{j}(2);
    end
end

% Find firing rate
rg=Range(XS,'s');
FR=length(Range(tsa))/(rg(end)-rg(1));

%% Plotting section

if plotresults
    
    figure('Color',[1 1 1])
    subplot(3,2,1), imagesc(map.time), axis xy, title('Occupation map')
    subplot(3,2,2), imagesc(map.count), axis xy, title('Spike map')
    subplot(3,2,3), imagesc(map.rate), axis xy, title('Firing map'), colorbar
    
    subplot(3,2,4) , plot(Data(XS),Data(YS),'Color',[0.8 0.8 0.8])
    hold on, plot(px,py,'r.')
    minx = min(Data(XS)); maxx = max(Data(XS)); difx = maxx-minx;
    xlim([minx-difx*0.05 maxx+difx*0.05]);
    miny = min(Data(YS)); maxy = max(Data(YS)); dify = maxy-miny;
    ylim([miny-dify*0.05 maxy+dify*0.05])
    title(['Firing rate : ',num2str(FR),' Hz'])
    
    subplot(3,2,5:6), hold on,
    plot(Range(XS,'s'),Data(XS))
    plot(Range(YS,'s'),Data(YS),'g')
    if maxy > maxx
        plot(Range(tsa,'s'),1.5*maxy*ones(length(Range(tsa,'s')),1),'r.')
    else
        plot(Range(tsa,'s'),1.3*maxx*ones(length(Range(tsa,'s')),1),'r.')
    end
    
end

if plotpoisson
    
    figure ('units', 'normalized','outerposition', [0 1 0.6 0.6])
    
    [A,B,C,D,E,G,F,H,I,J,map2,mapS2,stats2,px2,py2,FR2,sizeFinal2,PrField2,C2,ScField2, Ts]=...
        PlaceFieldPoisson(tsa, XS, YS, 'smoothing', smo, 'size', sizeMap, 'plotresults', 0);
    
    
    subplot(221)
    imagesc(map.rate), axis xy, title('Firing map'), colorbar
    c = clim;
    colormap jet
    subplot(222)
    plot(Data(XS),Data(YS),'Color',[0.8 0.8 0.8])
    hold on, plot(px,py,'r.')
    minx = min(Data(XS)); maxx = max(Data(XS)); difx = maxx-minx;
    xlim([minx-difx*0.05 maxx+difx*0.05]);
    miny = min(Data(YS)); maxy = max(Data(YS)); dify = maxy-miny;
    ylim([miny-dify*0.05 maxy+dify*0.05])
    rg=Range(XS,'s');
    title(['Firing rate : ',num2str(length(Range(tsa))/(rg(end)-rg(1))),' Hz'])
    subplot(223)
    imagesc(map2.rate), axis xy, title('Firing map Poisson'), colorbar
    colormap jet
    caxis([c]);
    subplot(224)
    plot(Data(XS),Data(YS),'Color',[0.8 0.8 0.8])
    hold on, plot(px2,py2,'r.')
    minx = min(Data(XS)); maxx = max(Data(XS)); difx = maxx-minx;
    xlim([minx-difx*0.05 maxx+difx*0.05]);
    miny = min(Data(YS)); maxy = max(Data(YS)); dify = maxy-miny;
    ylim([miny-dify*0.05 maxy+dify*0.05])
    title(['Firing rate : ',num2str(FR),' Hz'])
    annotation('textbox', [0.42 0.94 0.2 0.05], 'String', ['Spatial Info: ' num2str(stats.spatialInfo)], 'FontWeight', 'Bold',...
        'FitBoxToText','on', 'LineStyle','none', 'FontSize', 14)
    annotation('textbox', [0.42 0.48 0.2 0.05], 'String', ['Spatial Info: ' num2str(stats2.spatialInfo)], 'FontWeight', 'Bold',...
        'FitBoxToText','on', 'LineStyle','none', 'FontSize', 14)
end

%%
% ------------------------------- Helper functions -------------------------------


% CleanMap - repeatedly find and discard small regions of elevated firing rates

function map = CleanMap(M,N,threshold,map,n)

minFieldSize=5;

if nargin < 5
    n = 1;
elseif n >= 50
    % Maximum 50 recursive calls
    return
else
    n = n + 1;
end

peak = max(max(map));
peakLocation = FindPeakLocation(map);
field = FindFiringField(peakLocation(1),peakLocation(2),M,N,threshold*peak,map);
fieldSize = sum(sum(field));

if fieldSize > 0 && fieldSize < minFieldSize
    
    % Remove this region from rate map for next iteration
    map(logical(field)) = 0;
    map = CleanMap(M,N,threshold,map,n);
end

% FindPeakLocation - find the coordinates of the peak firing rate

function peakLocation = FindPeakLocation(map)

peak = max(max(map));
xy = (map == peak);
x = find(max(xy));
y = find(max(xy,[],2));
peakLocation = [x(1) y(1)];



