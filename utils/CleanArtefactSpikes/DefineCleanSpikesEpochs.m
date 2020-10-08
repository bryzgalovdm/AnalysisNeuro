function CleanEpoch = DefineCleanSpikesEpochs(file, Dir, StimDur)

%% Hyperparameters
Threshold = 1500;

%% Load data
load(file);
d = Range(raw);
TotalEpoch = intervalSet(d(1),d(end));

%% Check if you have information about stimulation
try
    load([Dir '/behavResources.mat'], 'StimEpoch');
    StimEpoch = intervalSet(Start(StimEpoch), Start(StimEpoch) + (StimDur+0.005)*1e4);
catch
    warning('No stim epoch. Will be thresholding only for removing atrifactual spikes');
end

%% Get high amplitude epoch
HIEpoch = thresholdIntervals(raw,Threshold,'Direction','Above');

%% Get final epoch
if exist('StimEpoch','var')
    CleanEpoch = TotalEpoch - or(StimEpoch, HIEpoch);
else
    CleanEpoch = TotalEpoch - HIEpoch;
end

end