%
% Asks if matrix M contains at least one value in vector "condition" at each position.
%
% USAGE:
%   result = iselement(M, condition)
%
% EXAMPLE: 
%   iselement(magic(5),1:10)
% 
%   ans =
%          0     0     1     1     0
%          0     1     1     0     0
%          1     1     0     0     0
%          1     0     0     0     1
%          0     0     0     1     1

%---------------------------- 
% Dirk Albrecht 
% Version 1.0 
% 22-May-2008 17:10:44 
%---------------------------- 

function result = iselement(M, condition)

valid = zeros(size(M));
for c = 1:length(condition)
    valid = valid | (M == condition(c));
end

result = valid;
