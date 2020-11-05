function ids = FindSessionID_ERC(behavResources, SessionName)
%
% Finds ids of a particular session types in ERC experiments
%
% INPUT
%
%     behavResources      structure of concatenated behavioral data
%     SessionName         Name to search
% 
%  OUTPUT
%
%     ids                 ids found
%
% Coded by Dima Bryzgalov, MOBS team, Paris, France
% 04/11/2020
% github.com/bryzgalovdm


ids = zeros(1,length(behavResources));

for k=1:length(behavResources)
    if ~isempty(strfind(behavResources(k).SessionName, SessionName))
        ids(k) = 1;
    end
end
ids=find(ids);

end