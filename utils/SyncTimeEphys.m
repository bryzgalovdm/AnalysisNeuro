function SyncTimeEphys(folderEvents, folderBehav, ExpeInfo)

% This funtion will synchronize the time recorded in Open Ephys with matlab-recorded
% behavior.
% You will need to provide the folder with Open-Ephys-recorded events (.mat)
% and ExpeInfo with correct info about digital inputs
% 
%  INPUT
% 
%     folder          folder with Open-Ephys-recorded events (.mat)
%     ExpeInfo        ExpeInfo with correct info about digital inputs
%
% Coded by Dima Bryzgalov, MOBS team, Paris, France
% 23/06/2020
% github.com/bryzgalovdm
    
%% Do the job
[~, out_ind] = regexp(folderEvents, 'recording\d/');
if ~exist([folderEvents(1:out_ind) 'continuous/continuous_Rhythm_FPGA-100.0.mat'], 'file')
    error('No syncronizing file. Please convert timestamps.npy in continuous folder');
end

TTLInfo_sess = MakeData_TTLInfo_OpenEphys([folderEvents 'events/Rhythm_FPGA-100.0_TTL_1.mat'],...
    folderEvents(1:out_ind),... 
    folderEvents(1:out_ind),...
    ExpeInfo);
% TTLInfo_sess = MakeData_TTLInfo_OpenEphys([folderEvents '/Rhythm_FPGA-100.0_TTL_1.mat'],...
%     [folderEvents(1:out_ind) 'continuous/Rhythm_FPGA-100.0/continuous_Rhythm_FPGA-100.0.mat'],...
%     ExpeInfo);

a = load([folderBehav '/behavResources.mat']);
a = CorrectTimeBackToEphys(a, TTLInfo_sess);

%% Save
mkdir([folderBehav '/old'])
if exist([folderBehav '/behavResources.mat'], 'file')
    movefile([folderBehav '/behavResources.mat'], [folderBehav '/old']);
end
save([folderBehav '/behavResources.mat'],'-struct','a')

end
