function idx = GetSessionIdx(behavResources, NameToSearch)
%
%  Finds indices of sessions the names of which contain a string NameToSearch
% 
% INPUT
% 
%     behavResources        structure containing all sessions
%     NAmeToSearch          a string to search into the structure field
%                           SessionName
% 
%  OUTPUT
%
%     idx                   indices in the structure
%
% by Dima Bryzgalov from MOBS team, Paris, France
% 29/05/2020
% github.com/bryzgalovdm

tempidx = zeros(length(behavResources),1);
for k=1:length(behavResources)
        if ~isempty(strfind(behavResources(k).SessionName, NameToSearch))
            tempidx(k) = 1;
        end
end
idx=find(tempidx);

% Issue the warning in case you did not find
if isempty(idx)
    warning(['No instances of ' NameToSearch ' were found']);
end

end