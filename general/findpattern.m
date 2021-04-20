%
% function index = findpattern(mat,testrow)
% 
%       finds rows in m-x-n matrix mat that match the criteria specified by either: 
%               1) a 1-x-n testrow vector specifying values exactly, or
%               2) a 2-x-n testrow matrix specifying the upper and lower
%               bounds for each column
%       use NaN to ignore a column and select all values for that column

%---------------------------- 
% Dirk Albrecht 
% Version 1.0 
% 16-Jun-2010 13:37:16 
%---------------------------- 

function index = findpattern(mat,testrow)

if size(testrow,1) == 1 testrow = [testrow;testrow]; end

match = [];
for col = 1:size(mat,2)
    if isnan(testrow(1,col))
        colmatch = ones(size(mat,1),1);
    else
        colmatch = mat(:,col)>=testrow(1,col) & mat(:,col)<=testrow(2,col);
    end
    match = [match, colmatch];
end

index = find(all(match,2));