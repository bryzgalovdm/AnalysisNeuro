function [Lfinal,LfSpk] = PlaceCellCrossValidation(Dir, cellnum, varargin)

%%% PlaceCellCrossValidation
%
% Estimating place cell spatial information by a cross-validating procedure
% Following Harris et al., Nature, 2003
%
% Dmitri Bryzgalov 2020 | Mobs team, France | github.com/bryzgalovdm
% ADAPTED FROM HD information
% Adrien Peyrache 2017 | peyrachelab.com | github.com/PeyracheLab
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

%% Parse optional parameters

for i=1:2:length(varargin)
    
    switch(lower(varargin{i}))
        
        case 'smoothing'
            smo = varargin{i+1};
            if ~isa(smo,'numeric')
                error('Incorrect value for property ''Smoothing'' (type ''help PlaceCellCrossValidation'' for details).');
            end
            
        case 'sizemap'
            sizeMap = varargin{i+1};
            if ~isa(sizeMap,'numeric')
                error('Incorrect value for property ''SizeMap'' (type ''help PlaceCellCrossValidation'' for details).');
            end
            
        case 'Verbose'
            verbose = varargin{i+1};
            if ~isa(sizeMap,'logical')
                error('Incorrect value for property ''Verbose'' (type ''help PlaceCellCrossValidation'' for details).');
            end
    end
end

% Default values of optional arguments

try
    smo;
catch
    smo=1;
end

try
    sizeMap;
catch
    sizeMap=50;
end

try
    verbose;
catch
    verbose=true;
end


%% Get data
load([Dir 'SpikeData.mat'], 'S');
load([Dir 'behavResources.mat'], 'SessionEpoch', 'CleanAlignedXtsd', 'CleanAlignedYtsd');

%% Treat the data

% Remove nans from position vars
% X
tempdatX = Data(CleanAlignedXtsd);
temptimeX = Range(CleanAlignedXtsd);
temptimeX(isnan(tempdatX)) = [];
tempdatX(isnan(tempdatX)) = [];
CleanAlignedXtsd = tsd(temptimeX, tempdatX);
% Y
tempdatY = Data(CleanAlignedYtsd);
temptimeY = Range(CleanAlignedYtsd);
temptimeY(isnan(tempdatY)) = [];
tempdatY(isnan(tempdatY)) = [];
CleanAlignedYtsd = tsd(temptimeY, tempdatY);

% Get time
Ts = End(SessionEpoch.Hab, 's') - Start(SessionEpoch.Hab, 's'); %%% in s
Td = End(SessionEpoch.Hab) - Start(SessionEpoch.Hab); %%% in tsd
time = Range(Restrict(CleanAlignedXtsd, SessionEpoch.Hab)); %%% in tsd
ddTs = Ts/length(time);


% We want to do a 10-fold cross-validation
% Length of 1/10th of the epoch
dT = Td/10;
dTs = Ts/10;

% Spike array
spk = Restrict(S{cellnum}, SessionEpoch.Hab);

%Likelihood function is divided in 10 values
Lf = zeros(10,1);

% Number of spikes
N = length(Data(spk));

% Average firing rate
fr = N/Ts;

%% Let's look at the information of a cell

for ep=1:10
    
    %Range of test set (offset by beginning of recording)
    tStart  = (ep-1)*dT + time(1);
    tEnd    = ep*dT + time(1);
    testIS = intervalSet(tStart, tEnd);
    
    % Defining training set (everything except test set)
    xTraining = Restrict(CleanAlignedXtsd, (SessionEpoch.Hab-testIS));
    yTraining = Restrict(CleanAlignedYtsd, (SessionEpoch.Hab-testIS));
    spkTraining = Restrict(spk, (SessionEpoch.Hab-testIS));
    
    % Estimating HD tuning curve on the training set
    [map, ~, ~, ~, ~, ~, xB, yB]=PlaceField_DB(spkTraining,xTraining, yTraining,...
        'SizeMap', sizeMap, 'Smoothing', smo, 'LargeMatrix', false, 'PlotResults', 0, 'PlotPoisson', 0);
    
    % Defining test set
    xTest = Restrict(CleanAlignedXtsd, testIS);
    yTest = Restrict(CleanAlignedYtsd, testIS);
    spkTest = Restrict(spk, testIS);
    
    % index of angBins corresponding to each angle during the test
    normX = (Data(xTest) - min(Data(xTest)))/(max(Data(xTest)) - min(Data(xTest)));
    normY = (Data(yTest) - min(Data(yTest)))/(max(Data(yTest)) - min(Data(yTest)));
    xx = floor((length(xB)-1)*normX)+1;
    yy = floor((length(yB)-1)*normY)+1;
    
    % expected firing rate (intensity function)
    intensityFct = diag(map.rate(xx,yy));
    
    % Evaluation of the intensity function at spike times
    pX = Restrict(xTest,spkTest,'align','closest');
    pY = Restrict(yTest,spkTest,'align','closest');
    
    dpX = Range(pX);
    dpY = Range(pY);
    idxX = zeros(length(pX),1);
    idxY = zeros(length(pY),1);
    for i=1:length(pX)
        idxX(i) = find(Range(xTest)==dpX(i));
        idxY(i) = find(Range(yTest)==dpY(i));
    end
    
    % Likelihood function
    % 'Flat' firing rate
    firingRate_flat = length(Data(spkTest))/dTs;
    
    % Log terms for the intensity functions
    logTermF = log2(intensityFct(idxX));
    logTermF(intensityFct(idxX)==0) = 0;
    
    % First term is very small; intensityFct is always <<< firingrate_flat. Why????
    Lf(ep) = - sum(intensityFct - firingRate_flat)*ddTs + sum(logTermF - log2(firingRate_flat));
    
end

Lfinal = sum(Lf)/Ts;
LfSpk = Lfinal/fr;

if verbose
    fprintf('Cross-validated information rate is %f (bit/sec)\n',Lfinal)
    fprintf('Cross-validated information per spike is %f (bit/spk)\n\n',LfSpk)
end

end