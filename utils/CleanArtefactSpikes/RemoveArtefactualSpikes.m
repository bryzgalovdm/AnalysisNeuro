function RemoveArtefactualSpikes(Dir, CleanEpoch)

%% Hyperparameters
samplingrate = 20000;
samplesinspike = 32;

%% Create folder to store original files
if exist([Dir '/orig_spikefiles'], 'dir') ~= 7
    mkdir([Dir '/orig_spikefiles']);
end
    

%% Load ExpeInfo
load([Dir '/ExpeInfo.mat']);
BaseFileName = ['M' num2str(ExpeInfo.nmouse) '_' ExpeInfo.date '_' ExpeInfo.SessionType];


for igroup = 1:ExpeInfo.SpikeGroupInfo.SpikeGroupNum
    %% Load the data
    ChanNum = length(str2num(ExpeInfo.SpikeGroupInfo.ChanNames{igroup}));
    
    resfile = importdata([Dir '/' BaseFileName '.res.' num2str(igroup)]);
    
    spkfid = fopen([Dir '/' BaseFileName '.spk.' num2str(igroup)]);
    spkfile = fread(spkfid, 'int16');
    spkfile = reshape(spkfile, [ChanNum, samplesinspike,length(spkfile)/(ChanNum*samplesinspike)]);

    %% Remove spikes with artefacts
    
    rests = ts(resfile);
    CleanEpoch_samples = intervalSet(Start(CleanEpoch)/1e4*samplingrate, End(CleanEpoch)/1e4*samplingrate);
    
    [newresfile, idxs] = Restrict(rests, CleanEpoch_samples);
    newresfile = Range(newresfile);
    
    newspkfile = spkfile(:,:,idxs);
    newspkfile = newspkfile(:);
    
    %% Move original ones to <old> folder and save new ones
    
    % Move <old> ones
    movefile([Dir '/' BaseFileName '.res.' num2str(igroup)], [Dir '/orig_spikefiles']);
    movefile([Dir '/' BaseFileName '.spk.' num2str(igroup)], [Dir '/orig_spikefiles']);
    
    % Save new ones
    dlmwrite([Dir '/' BaseFileName '.res.' num2str(igroup)],...
        newresfile, 'delimiter', '\n', 'precision', 40);
    
    fspk = fopen([Dir '/' BaseFileName '.spk.' num2str(igroup)], 'w');
    fwrite(fspk, newspkfile, 'int16');
    fclose(fspk);
    

end