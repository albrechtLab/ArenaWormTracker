%
% Split input matrix into n parts along dim1, adding split sections along dim2.
%
% USAGE:
%     output = split(input,n,dim1,dim2)
%
% EXAMPLE:
%     a = repmat((1:3)' * (1:3),[1 1 2])
% 
%     a(:,:,1) =
%          1     2     3
%          2     4     6
%          3     6     9
% 
% 
%     a(:,:,2) =
%          1     2     3
%          2     4     6
%          3     6     9
% 
%     split(a,2,3,1)
% 
%     ans =
%          1     2     3
%          2     4     6
%          3     6     9
%          1     2     3
%          2     4     6
%          3     6     9

%---------------------------- 
% Dirk Albrecht 
% Version 1.0 
% 14-Apr-2010 22:51:05 
%---------------------------- 

function output = split(input,n,dim1,dim2)

sizes = size(input);
delta = ceil(sizes(dim1)/n);
maxdims = max([ ndims(input), dim1, dim2 ]);

dimshift = 1:maxdims; dimshift(dim1)=1; dimshift(1)=dim1;
input = permute(input,dimshift);

output = [];
for i = 1:n
    index = (i-1)*delta+1:i*delta;
    output = safecat(dimshift(dim2),output,input(index,:,:,:));
end

dimshift = 1:ndims(output); dimshift(dim1)=1; dimshift(1)=dim1;
output = permute(output,dimshift);