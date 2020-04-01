function [map,mapS,stats,px,py,FR,sizeFinal,PrField,C,ScField,pfH,pf]=PlaceField_DB(tsa,XS,YS,varargin)



%  OUTPUT
%
%    map.x          x bins
%    map.y          y bins
%    map.rate       average firing rate map (in Hz)
%    map.time       occupancy map (in s)
%    map.count      firing count map
%    stats.x        abscissa of the peak (in bins)
%    stats.y        ordinate of the peak (in bins)
%    stats.peak     in-field peak firing rate (in Hz)
%    stats.mean     in-field mean firing rate (in Hz)
%    stats.size     firing field size (in bins)
%    stats.field    firing field (1 = bin in field, 0 = bin not in field)
%    stats.fieldX   firing field x boundaries (in bins)
%    stats.fieldY   firing field y boundaries (in bins)
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

oldversion=0;
LMatrix=1;


for i=1:2:length(varargin)
    
    %           if ~isa(varargin{i},'char'),
    %            error(['Parameter ' num2str(i) ' is not a property (type ''help ICSSexplo'' for details).']);
    %         end
    
    switch(lower(varargin{i})),
        
        case 'smoothing',
            smo = varargin{i+1};
            if ~isa(smo,'numeric'),
                error('Incorrect value for property ''smoothing'' (type ''help PlaceField'' for details).');
            end
            
        case 'figure',
            plo = varargin{i+1};
            if ~isa(plo,'numeric'),
                error('Incorrect value for property ''figure'' (type ''help PlaceField_DB'' for details).');
            end
        case 'video',
            freqVideo = varargin{i+1};
            if ~isa(freqVideo,'numeric'),
                error('Incorrect value for property ''video'' (type ''help PlaceField_DB'' for details).');
            end
            
        case 'threshold',
            threshold = varargin{i+1};
            if ~isa(threshold,'numeric'),
                error('Incorrect value for property ''threshold'' (type ''help PlaceField_DB'' for details).');
            end
            
        case 'size',
            sizeMap = varargin{i+1};
            if ~isa(sizeMap,'numeric'),
                error('Incorrect value for property ''size'' (type ''help PlaceField_DB'' for details).');
            end
            
        case 'limitmaze',
            lim = varargin{i+1};
            
        case 'largematrix',
            LMatrix = varargin{i+1};
            
        case 'plotresults',
            plotresults = varargin{i+1};
            if ~(plotresults == 1 || plotresults == 0)
                error('Incorrect value for property ''plotresults'' (type ''help PlaceField_DB'' for details).');
            end
            
        case 'plotpoisson',
            plotpoisson = varargin{i+1};
            if ~(plotpoisson == 1 || plotpoisson == 0)
                error('Incorrect value for property ''plotpoisson'' (type ''help PlaceField_DB'' for details).');
            end
            
    end
end




try
    threshold;
catch
    threshold = 0.7;
end


try
    plo;
catch
    plo=0;
end



try
    smo;
catch
    smo=3;
end


try
    sizeMap;
    
catch
    sizeMap=50;
    
end

try
    freqVideo;
catch
    freqVideo=30;
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


px =Data(Restrict(XS,tsa,'align','closest'));
py =Data(Restrict(YS,tsa,'align','closest'));


if 0
    figure('Color',[1 1 1]), plot(Data(XS),Data(YS),'Color',[0.8 0.8 0.8])
    hold on, plot(px,py,'r.')
    xlim([min(Data(XS)) max(Data(XS))])
    ylim([min(Data(YS)) max(Data(YS))])
end

try
    lim(1:2)=sort(lim(1:2));
    lim(3:4)=sort(lim(3:4));
end



if oldversion
    
    try
        try
            [occH, x1, x2] = hist2d([lim(1) ;Data(XS) ;lim(2)], [lim(3) ;Data(YS) ;lim(4)], sizeMap, sizeMap);
        catch
            [occH, x1, x2] = hist2d([lim(1) ;Data(XS) ;lim(2)], [lim(1) ;Data(YS) ;lim(2)], sizeMap, sizeMap);
        end
        occH(1,1)=0;
        occH(end,end)=0;
    catch
        [occH, x1, x2] = hist2d(Data(XS), Data(YS), sizeMap, sizeMap);
    end
    
else
    
    
    try
        try
            pas12=(lim(2)-lim(1))/(sizeMap-2);
            pas34=(lim(4)-lim(3))/(sizeMap-2);
            
            [occH, x1, x2] = hist2d(Data(XS), Data(YS), lim(1)-pas12/2:pas12:lim(2)+pas12/2, lim(3)-pas34/2:pas34:lim(4)+pas34/2);
        catch
            pas12=(lim(2)-lim(1))/(sizeMap-2);
            [occH, x1, x2] = hist2d(Data(XS), Data(YS), lim(1)-pas/2:pas:lim(2)+pas/2, lim(1)-pas12/2:pas12:lim(2)+pas12/2);
        end
        
    catch
        [occH, x1, x2] = hist2d(Data(XS), Data(YS), sizeMap, sizeMap);
    end
    
    
    
end



sz = length(occH(:));
sOcc = sort(occH(:));
th = sOcc(end-(ceil(sz/100)));

pX = Restrict(XS,tsa,'align','closest');
pY = Restrict(YS, tsa,'align','closest');

sigmaS = 10;



if oldversion
    
    
    try
        pfH = hist2d([lim(1) ;Data(pX) ;lim(2)], [lim(1) ;Data(py) ;lim(2)], x1, x2);
        pfH(1,1)=0;
        pfH(end,end)=0;
    catch
        pfH = hist2d(Data(pX), Data(pY), x1, x2);
    end
    
    
else
    
    try
        try
            pas12=(lim(2)-lim(1))/(sizeMap-2);
            pas34=(lim(4)-lim(3))/(sizeMap-2);
            pfH = hist2d(Data(pX), Data(py), lim(1)-pas12/2:pas12:lim(2)+pas12/2, lim(3)-pas34/2:pas34:lim(4)+pas34/2);
        catch
            
            pas12=(lim(2)-lim(1))/(sizeMap-2);
            pfH = hist2d(Data(pX), Data(py), lim(1)-pas/2:pas:lim(2)+pas/2, lim(1)-pas12/2:pas12:lim(2)+pas12/2);
            
            
        end
        
        
    catch
        pfH = hist2d(Data(pX), Data(pY), x1, x2);
    end
    
    
end








%     dbclear warning
%     warning off

pf = 30 * pfH./occH;

if LMatrix
    largerMatrix = zeros(sizeMap+floor(sizeMap/4),sizeMap+floor(sizeMap/4));
    largerMatrix(1+floor(sizeMap/8):sizeMap+floor(sizeMap/8),1+floor(sizeMap/8):sizeMap+floor(sizeMap/8)) = pf;
    % largerMatrix = zeros(250,250);
    % largerMatrix(26:225,26:225) = pf;
    pf = largerMatrix;
    
end


warning on
pf(isnan(pf))=0;
%
% for x=1:smo;
% for y=1:smo;
%     gw(x,y) = exp(-((x-5)*(x-5) + (y-5)*(y-5))/(2*sigmaS));
% end;
% end

sg = sort(pf(~isnan(pf(:))));
th = sg(end-5);
pf(pf>th) = 0;

%    pf = conv2(pf,gw,'same')/sum(sum(gw));
pf=SmoothDec(pf,[smo,smo]);



if 0
    figure, clf
    
    fh = imagesc(pf');
    colormap(jet)
    axis xy
    colorbar
end

pf=pf';
map.rate=pf;

%%%%%%%%%%%%%%%occH = conv2(occH,gw,'same')/sum(sum(gw));

%occH=SmoothDec(occH,[smo,smo]);
if LMatrix
    
    largerMatrix = zeros(sizeMap+floor(sizeMap/4),sizeMap+floor(sizeMap/4));
    largerMatrix(1+floor(sizeMap/8):sizeMap+floor(sizeMap/8),1+floor(sizeMap/8):sizeMap+floor(sizeMap/8)) = occH';
    largerMatrix=SmoothDec(largerMatrix,[smo,smo]);
    
    % largerMatrix = zeros(250,250);
    % largerMatrix(26:225,26:225) = occH';
    map.time=largerMatrix/freqVideo;
else
    
    occH=SmoothDec(occH'/freqVideo,[smo,smo]);
    map.time=occH;
end



%%%%%%%%%%%%%%%%pfH = conv2(pfH,gw,'same')/sum(sum(gw));
%pfH=SmoothDec(pfH,[smo,smo]);
if LMatrix
    
    largerMatrix = zeros(sizeMap+floor(sizeMap/4),sizeMap+floor(sizeMap/4));
    largerMatrix(1+floor(sizeMap/8):sizeMap+floor(sizeMap/8),1+floor(sizeMap/8):sizeMap+floor(sizeMap/8)) = pfH';
    largerMatrix=SmoothDec(largerMatrix,[smo,smo]);
    
    % largerMatrix = zeros(250,250);
    % largerMatrix(26:225,26:225) = pfH';
    map.count=largerMatrix;
    
else
    
    pfH=SmoothDec(pfH',[smo,smo]);
    map.count=pfH;
end

%--------------------------------------------------------------
%--------------------------------------------------------------
%--------------------------------------------------------------

%mapS.time=SmoothDec(map.time,[smo/2 smo/2]);
%mapS.count=SmoothDec(map.count,[smo/2 smo/2]);
%mapS.rate=SmoothDec(map.rate,[smo/2 smo/2]);

%--------------------------------------------------------------
%--------------------------------------------------------------
%--------------------------------------------------------------

%threshold = 0.7;

% nBins = 250;

if LMatrix
    nBins =sizeMap+floor(sizeMap/4);
else
    nBins =sizeMap;
end


nBinsX = nBins(1);
if length(nBins) == 1,
    nBinsY = nBinsX;
    nBins(2) = nBins;
else
    nBinsY = nBins(2);
end

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
if T == 0,
    stats.specificity = 0;
else
    occupancy = map.time/(T+eps);
    meanFiringRate = sum(sum(map.count))/(sum(sum(map.time)+eps));
    if meanFiringRate == 0,
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

if max(max(map.rate)) == 0,
    stats.field = logical(zeros(size(map.rate)));
    return;
end

%keyboard
% First clean the map from very small regions of elevated firing rates
map.rate = CleanMap(nBinsX,nBinsY,threshold,map.rate,20);

mapS=map;
pb=mapS.time;
fr=mapS.rate;
fr(pb<sum(sum(pb))/10000)=0;
mapS.rate=fr;

% In the absence of any information regarding the location of the center,
% find two candidate firing fields
rate = map.rate;
for i = 1:2,
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
if fieldSize(2) == 0,
    winner = 1;
else
    %% peakRelativeDifference = abs(peak(2)-peak(1))/((peak(1)+peak(2))/2);
    % 		sizeRelativeDifference = (fieldSize(2)-fieldSize(1))/((fieldSize(1)+fieldSize(2))/2);
    sizeRelativeDifference = fieldSize(2)/fieldSize(1);
    if sizeRelativeDifference > 3, % Used to be 0.2
        winner = 2; % Choose firing field #2 (see below)
    elseif sizeRelativeDifference < 0.33, % used to be else
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
    if ~isempty(i),
        stats.fieldX = [i(1) i(end)];
    end
    i = find(max(stats.field,[],2));
    if ~isempty(i),
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
        if ~isempty(j),
            stats.fieldX{i} = [j(1) j(end)];
        end
        j = find(max(stats.field{i},[],2));
        if ~isempty(j),
            stats.fieldY{i} = [j(1) j(end)];
        end
    end
end

rg=Range(XS,'s');
idxx=(find(diff(Range(XS,'ms'))>10*median(diff(Range(XS,'ms')))));


if idxx>0
    for i=1:length(idxx)
        %timeEpoch(i)=rg(end)-rg(idxx(i)+1)+ rg(idxx(i)-1)-rg(1);
        timeEpoch(i)=rg(end)-rg(idxx(i)+1)+ rg(max(2,idxx(i))-1)-rg(1); % add by Marie on 16January2014
    end
    
    timeEpoch=sum(timeEpoch);
else
    timeEpoch=rg(end)-rg(1);
end


meanrate=length(tsa)/timeEpoch;

%stat=SpatialInfo(map,meanrate);
stat=SpatialInfo(map);

stats.spatialInfo=stat.info;
stats.sparsity=stat.sparsity;

MT=map.time;
MR=map.rate;
iFr=MR(find(MT>1E-2));

Mt=zeros(size(MT,1),size(MT,2));
for l=1:size(MT,1)
    for j=1:size(MT,2)
        if MT(l,j)>1E-2
            Mt(l,j)=1;
        end
    end
end

if iscell(stats.field)
    for j=1:length(stats.field)
        C{j}=GravityCenter(stats.field{j}.*map.rate);
        GC{j}(j,1)=C{j}(1);
        GC{j}(j,2)=C{j}(2);
    end
end

if plo==1
    if iscell(stats.field)
        if length(stats.field)>1
            fieldRR = stats.field{1}+stats.field{2};
            figure('Color',[1 1 1]), imagesc(fieldRR.*map.rate+Mt), axis xy
            hold on, plot(C(1),C(2),'wo','MarkerFaceColor','k')
            title(['Firing Map : Field size: ',num2str(stats.size)])
        else
            figure('Color',[1 1 1]), imagesc(stats.field.*map.rate+Mt), axis xy
            hold on, plot(C(1),C(2),'wo','MarkerFaceColor','k')
            title(['Firing Map : Field size: ',num2str(stats.size)])
        end
    end
end

if plotresults
    
    figure('Color',[1 1 1]),
    %clf, subplot(3,2,1), imagesc(SmoothDec(map.time,[sizeMap/50,sizeMap/50])), axis xy, title('Occupation map')
    %subplot(3,2,2), imagesc(SmoothDec(map.count,[sizeMap/50,sizeMap/50])), axis xy, title('Spike map')
    %subplot(3,2,3), imagesc(SmoothDec(map.rate,[sizeMap/50,sizeMap/50])), axis xy, title('Firing map'), colorbar
    
    
    % clf, subplot(3,2,1), imagesc(SmoothDec(map.time,[smo smo])), axis xy, title('Occupation map')
    % subplot(3,2,2), imagesc(SmoothDec(map.count,[smo smo])), axis xy, title('Spike map')
    % subplot(3,2,3), imagesc(SmoothDec(map.rate,[smo smo])), axis xy, title('Firing map'), colorbar
    
    clf, subplot(3,2,1), imagesc(map.time), axis xy, title('Occupation map')
    subplot(3,2,2), imagesc(map.count), axis xy, title('Spike map')
    subplot(3,2,3), imagesc(map.rate), axis xy, title('Firing map'), colorbar
    
    
    
    subplot(3,2,4) , plot(Data(XS),Data(YS),'Color',[0.8 0.8 0.8])
    hold on, plot(px,py,'r.')
    xlim([min(Data(XS))-5 max(Data(XS))+5])
    ylim([min(Data(YS))-5 max(Data(YS))+5])
    rg=Range(XS,'s');
    title(['Firing rate : ',num2str(length(Range(tsa))/(rg(end)-rg(1))),' Hz'])
    
    subplot(3,2,5:6), hold on,
    plot(Range(XS,'s'),Data(XS))
    plot(Range(YS,'s'),Data(YS),'g')
    plot(Range(tsa,'s'),100*ones(length(Range(tsa,'s')),1),'r.')
    
end

FR=length(Range(tsa))/(rg(end)-rg(1));

PrField=field{1};
ScField=field{2};

if LMatrix
    sizeFinal=sizeMap+floor(sizeMap/4);
    
else
    
    sizeFinal=sizeMap;
    
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
    xlim([min(Data(XS))-0.1 max(Data(XS))+0.1])
    ylim([min(Data(YS))-0.1 max(Data(YS))+0.1])
    rg=Range(XS,'s');
    title(['Firing rate : ',num2str(length(Range(tsa))/(rg(end)-rg(1))),' Hz'])
    subplot(223)
    imagesc(map2.rate), axis xy, title('Firing map Poisson'), colorbar
    colormap jet
    caxis([c]);
    subplot(224)
    plot(Data(XS),Data(YS),'Color',[0.8 0.8 0.8])
    hold on, plot(px2,py2,'r.')
    xlim([min(Data(XS))-0.1 max(Data(XS))+0.1])
    ylim([min(Data(YS))-0.1 max(Data(YS))+0.1])
    rg=Range(XS,'s');
    title(['Firing rate : ',num2str(length(Range(Ts))/(rg(end)-rg(1))),' Hz'])
    annotation('textbox', [0.42 0.94 0.2 0.05], 'String', ['Spatial Info: ' num2str(stats.spatialInfo)], 'FontWeight', 'Bold',...
        'FitBoxToText','on', 'LineStyle','none', 'FontSize', 14)
    annotation('textbox', [0.42 0.48 0.2 0.05], 'String', ['Spatial Info: ' num2str(stats2.spatialInfo)], 'FontWeight', 'Bold',...
        'FitBoxToText','on', 'LineStyle','none', 'FontSize', 14)
end


% ------------------------------- Helper functions -------------------------------


% CleanMap - repeatedly find and discard small regions of elevated firing rates

function map = CleanMap(M,N,threshold,map,n)

minFieldSize=5;

if nargin < 5,
    n = 1;
elseif n >= 50,
    % Maximum 50 recursive calls
    return
else
    n = n + 1;
end

peak = max(max(map));
peakLocation = FindPeakLocation(map);
field = FindFiringField(peakLocation(1),peakLocation(2),M,N,threshold*peak,map);
fieldSize = sum(sum(field));
% disp(['CleanMap: field size = ' int2str(fieldSize)]);

% if fieldSize > 0 & fieldSize < SETTINGS.minFieldSize,

if fieldSize > 0 & fieldSize < minFieldSize,
    
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



