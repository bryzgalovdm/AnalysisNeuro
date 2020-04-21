function [Lfinal,LfSpk] = PlaceCellCrossValidation(ClusterTSD, Xtsd, Ytsd, varargin)

%%% PlaceCellCrossValidation
%
% INPUT
% 
%     ClusterTSD            spike times in tsd
%     Xtsd                  X position in tsd
%     Ytsd                  Y position is tsd
%     
%     optional arguments:
%     
%     Smoothing             spatial smooting factor (default - 0)
%     SizeMap               size of the map (it will be size * size) in pixels
%     Verbose               verbose result (default - true)


%  OUTPUT
%
%     Lfinal                Predictability in bits/sec
%     LfSpk                 Predictability in bits/spike
%     
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
            
        case 'verbose'
            verbose = varargin{i+1};
            if ~isa(verbose,'logical')
                error('Incorrect value for property ''Verbose'' (type ''help PlaceCellCrossValidation'' for details).');
            end
    end
end

% Default values of optional arguments

try
    smo;
catch
    smo=0;
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


%% Treat the data

% Remove nans from position vars
% X
tempdatX = Data(Xtsd);
temptimeX = Range(Xtsd);
temptimeX(isnan(tempdatX)) = [];
tempdatX(isnan(tempdatX)) = [];
CleanXtsd = tsd(temptimeX, tempdatX);
% Y
tempdatY = Data(Ytsd);
temptimeY = Range(Ytsd);
temptimeY(isnan(tempdatY)) = [];
tempdatY(isnan(tempdatY)) = [];
CleanYtsd = tsd(temptimeY, tempdatY);

% Get time
time = Range(Xtsd); %%% in tsd
Ts = (time(end) - time(1))/1e4; %%% in s
Td = time(end) - time(1); %%% in tsd
ddTs = Ts/length(time);
EpochFull = intervalSet(time(1),time(end));


% We want to do a 10-fold cross-validation
% Length of 1/10th of the epoch
dT = Td/10;
dTs = Ts/10;

% Spike array
spk = ClusterTSD;

%Likelihood function is divided in 10 values
Lf = zeros(10,1);

% Number of spikes
N = length(Data(spk));

% Average firing rate
fr = N/Ts;

% Spatial bins with FR <= <exlThresh> become 0 Hz firing rate in the log-term 
exlThresh = 1/dTs; % Hz

%% Let's look at the information of a cell

for ep=1:10
    
    %Range of test set (offset by beginning of recording)
    tStart  = (ep-1)*dT + time(1);
    tEnd    = ep*dT + time(1);
    testIS = intervalSet(tStart, tEnd);
    
    % Defining training set (everything except test set)
    xTraining = Restrict(CleanXtsd, (EpochFull-testIS));
    yTraining = Restrict(CleanYtsd, (EpochFull-testIS));
    spkTraining = Restrict(spk, (EpochFull-testIS));
    
    % Estimating place cell rate map on the training set
    [map, mapNS, ~, ~, ~, ~, xB, yB]=PlaceField_DB(spkTraining,xTraining, yTraining,...
        'SizeMap', sizeMap, 'LargeMatrix', false, 'PlotResults', 0, 'PlotPoisson', 0);
    
    % Defining test set
    xTest = Restrict(CleanXtsd, testIS);
    yTest = Restrict(CleanYtsd, testIS);
    spkTest = Restrict(spk, testIS);
    
    % index of X and Y coordinates corresponding to each spatial bin during the test
    normX = (Data(xTest) - min(Data(xTest)))/(max(Data(xTest)) - min(Data(xTest)));
    normY = (Data(yTest) - min(Data(yTest)))/(max(Data(yTest)) - min(Data(yTest)));
    xx = floor((length(xB)-1)*normX)+1;
    yy = floor((length(yB)-1)*normY)+1;
    
    % expected firing rate (intensity function)
    intensityFct = diag(map.rate(xx,yy));
    
    % Evaluation of the intensity function at spike times
    pX = Restrict(xTest,spkTest,'align','closest');
    
    dpX = Range(pX);
    idxX = zeros(length(pX),1);
    for i=1:length(pX)
        idxX(i) = find(Range(xTest)==dpX(i));
    end
    
    % Likelihood function
    % 'Flat' firing rate
    firingRate_flat = length(Data(spkTest))/dTs;
    
    % Log terms for the intensity functions
    logTermF = log2(intensityFct(idxX));
    logTermF(intensityFct(idxX)<exlThresh) = 0;
    
    % Do a final calculation
    Lf(ep) = - sum(intensityFct - firingRate_flat)*ddTs + sum(logTermF - log2(firingRate_flat));
    
end

Lfinal = sum(Lf)/Ts;
LfSpk = Lfinal/fr;

if verbose
    fprintf('Cross-validated information rate is %f (bit/sec)\n',Lfinal)
    fprintf('Cross-validated information per spike is %f (bit/spk)\n\n',LfSpk)
end

end