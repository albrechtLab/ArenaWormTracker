%
% Compare vectors
%
% USAGE:
%   [duplicate, aunique, bunique] = veccomp(a, b)
%   
%   a and b should be row vectors
%   zeros are not compared
%   duplicate returns NaN if no overlap

%---------------------------- 
% Dirk Albrecht 
% Version 1.0 
% 11-Feb-2008 12:37:38 
%---------------------------- 

function [duplicate, aunique, bunique] = veccomp(a, b)

test = zeros(size(b));
for i = 1:length(a)
    test = test + (b == a(i));
end

duplicate = b(find(test));
bunique = b(find(~test));
bunique = bunique(find(bunique ~= 0));

test = zeros(size(a));
for i = 1:length(b)
    test = test + (a == b(i));
end

%duplicate = a(find(test));
aunique = a(find(~test));
aunique = aunique(find(aunique ~= 0));

if length(duplicate) == 0 duplicate = NaN; end
if length(aunique) == 0 aunique = NaN; end
if length(bunique) == 0 bunique = NaN; end