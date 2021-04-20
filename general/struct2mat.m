%
% Extract data from structures.
%
% USAGE:
%   output = struct2mat(dim,structure,index,fields)
%   
%   dim: dimension along which to concatenate output matrix
%   structure: input struct
%   index: indices of structure to concatenate ([] for all, default)
%   fields: cell array of field names
%       for structure.field, use {'field'}
%       for structure.substruct.field, use {'substruct','field'}

%---------------------------- 
% Dirk Albrecht 
% Version 1.0 
% 02-Apr-2010 15:35:17 
%---------------------------- 

function output = struct2mat(dim,structure,index,fields)

if ~iscell(fields) fields = {fields}; end
if isempty(index) index = 1:numel(structure); end

output = [];
for i = 1:length(index)
    substructure = structure(index(i));
    for j = 1:length(fields)
        substructure = getfield(substructure,char(fields(j)));
    end
    
    output = safecat(dim,output,substructure);
end