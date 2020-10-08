function [stimtimes, stim_idx] = SimulateStim_PositionDecoder(DecodedData, EL_thresh, T_thresh)

%% Prepare arrays

% Apply threshold
data_idx = find(DecodedData(:,4)<EL_thresh);
data = DecodedData(data_idx,:);

%% Find stim indices
for it = 2:length(data)
    if data(it,2) > 0.8 && data(it,2) < 1 && data(it,3) > 0.42 && data(it,3) < 0.65
        if data(it-1,2) == data(it,2) &&  data(it-1,3) == data(it,3)
        temp(it-1) = it-1;
        end
    end
end
if exist('temp', 'var')
    stimIDX = nonzeros(temp);
    clear temp
else
   stimtimes = [];
   stim_idx = [];
   return
end

%% Allocate data
stimtimes = nan(length(stimIDX),1);
stim_idx = nan(length(stimIDX),1);

%% Get stimulation times and indices
tempdiff = diff(data(:,1));
stimtimes(1) = DecodedData(data_idx(stimIDX(1)),1); % Stimulation to the first stimulus
stim_idx(1) = data_idx(stimIDX(1)); % Stimulation to the first stimulus

count = 0;
for numstim = 2:length(stimIDX)
    count = count + sum(tempdiff(stimIDX(numstim-1):stimIDX(numstim)));
    if count > T_thresh
        count = 0;
        stimtimes(numstim) = DecodedData(data_idx(stimIDX(numstim)),1);
        stim_idx(numstim) = data_idx(stimIDX(numstim));
    end
end
clear tempdiff count

%% Clean resulting arrays
stimtimes(isnan(stimtimes)) = [];
stim_idx(isnan(stim_idx)) = [];

end