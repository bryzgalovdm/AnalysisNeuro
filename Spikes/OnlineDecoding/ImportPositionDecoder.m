function ImportPositionDecoder(filepath, oebin_folder, savefolder)

%% Learn sampling rate and start time
% Sampling rate
oebin = fileread([oebin_folder 'structure.oebin']);
[~, sr_id] = regexp(oebin,'"sample_rate": ');
samplingrate = str2double(oebin(sr_id(1)+1:sr_id(1)+5));

% Start time
sync = load([oebin_folder 'continuous/continuous_Rhythm_FPGA-100.0.mat']);
starttime = sync.timestamps(1);

% Legacy start time
% sync = fileread([sync_folder 'sync_messages.txt']);
% [~,sync_id_st] = regexp(sync,'start time: ');
% sync_id_en = regexp(sync,'@');
% sync_id_en = sync_id_en(2)-1;
% starttime = str2double(sync(sync_id_st:sync_id_en));

%% Load the file
a = load(filepath);

%% Sort
[timestamps, idx] = sort(a.timestamps);
% Apply sorting
channel_states = a.channel_states(idx);
channels = a.channels(idx);
full_words = a.full_words(idx);
metadata = a.metadata(idx,:);

%% Correct the time
time = double(timestamps - starttime)/samplingrate*1e4;
% Remove negative times
idx_neg = find(time<0);
time(idx_neg) = [];
channel_states(idx_neg) = [];
channels(idx_neg) = [];
full_words(idx_neg) = [];
metadata(idx_neg,:) = [];

stimtimes = time(channel_states==1);
DecodedStimIDX = find(channel_states==1);

%% Construct variables to save
DecodedXtsd = tsd(time, metadata(:,1));
DecodedYtsd = tsd(time, metadata(:,2));
DecodedELtsd = tsd(time, metadata(:,3));
DecodedData = [time'/1e4 metadata];

DecodedStimEpoch = intervalSet(stimtimes, stimtimes);

%% Save
save([savefolder '/behavResources.mat'], 'DecodedXtsd', 'DecodedYtsd', 'DecodedELtsd', 'DecodedData', 'DecodedStimEpoch',...
    'DecodedStimIDX', '-append');

end