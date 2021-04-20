%---------------------------- 
% Dirk Albrecht 
% Version 1.0 
% 01-Mar-2010 12:15:17 
%---------------------------- 

function output = scalenoisy(input,inputrange,outputrange)
if diff(size(input)) > 0 input = input'; end

if nargin < 3 outputrange = [0 1]; end
if nargin < 2 || isempty(inputrange) inputrange = 1:length(input); end
if max(inputrange) <= 1 inputrange = max(1,[round(inputrange(1)*length(input)):inputrange(2)*length(input)]); end

N = 10;
[a,b]=hist(input(inputrange,:),N^2); if ~any(size(a) == 1) a = sum(a')'; end
%divide into N bins
binmax = max(reshape(a,N,N));
lowbin = find(diff(binmax) <= -10,1);
highbin = N+1-find(diff(fliplr(binmax)) <= -10,1);

lowr = N*lowbin+(-N+1:0);
highr = N*highbin+(-N+1:0);

fullrange = [mean(b(min(lowr)-1+find(a(lowr) == max(a(lowr)),1))), ...
             mean(b(min(highr)-1+find(a(highr) == max(a(highr)),1)))]; 
z =([min(b) max(b)]-fullrange(1))/diff(fullrange);

output = (scale(input, z(1),z(2))) * diff(outputrange) + outputrange(1);