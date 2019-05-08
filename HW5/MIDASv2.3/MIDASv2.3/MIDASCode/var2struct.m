function out=var2struct(varargin)
for i=1:length(varargin)
    if isempty(inputname(i))
        error('You have to provide variables, not numbers or results of the functions.')
    end
    out.(inputname(i))=varargin{i};
end
end