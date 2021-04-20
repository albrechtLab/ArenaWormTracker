%
% Converts all numeric elements in a structure to single precision.
%
% USAGE:
%   output = singleStruct(input)

%---------------------------- 
% Dirk Albrecht 
% Version 1.0 
% 28-Mar-2011 18:45:26 
%---------------------------- 

function output = singleStruct(input)

fields = fieldnames(input);

output = input;
for i = 1:numel(input)
    for f = 1:length(fields)
        value = getfield(input,{i},fields{f});
        if isnumeric(value)
            output = setfield(output,{i},fields{f},single(value));
        else
            output = setfield(output,{i},fields{f},value);
        end
    end
end
