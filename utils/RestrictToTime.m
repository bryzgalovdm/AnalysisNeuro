function I = RestrictToTime(is, len)

% I = RestrictToTime(is, len) restricts IntervalSet to be of specific length
%
% This code makes your interval to be specific length
% 
% 
% USAGE
%
%    I = RestrictToTime(is, len)
%
%    INPUTS
%
%    is         an interval set
%    len        desired total length of the resulting interval set in t
%
%    OUTPUTS
%
%    I      resulting resricted interval set
%
% 
%    See also Restrict, intervalSet
%
% the current program is an addition to Francesco P. Battaglia's tsd plugin
% copyright (c) 2004 Francesco P. Battaglia, (c) 2018-2019 Dmitri Bryzgalov
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% released under the GPL   

% Find length of each interval
for i=1:length(Start(is))
    l=(End(is)-Start(is));
end

% Find the interval that exceeds the desired length
count=0;
j=0;
if sum(l)>len
    while count<len
        j=j+1;
        count = count + l(j);
    end
else
    error('Length of interval set is shorter than time you are trying to restrict it to');
end
last = j;

% Crop the interval that exceeds the desired length
if count > len
    almost = sum(l(1:last-1));
    delta = len - almost;
    
    is_almost = subset(is, 1:last-1);
    is_delta = intervalSet(Start(subset(is, last)), Start(subset(is, last))+delta);
    
    is_result = or(is_almost, is_delta);
elseif count == len
    is_result = subset(is, 1:last);
else
    error('The author is a fool. Let him know - dmitri.bryzgalov@espci.fr');
end

% Final assignment
I = is_result;

end


