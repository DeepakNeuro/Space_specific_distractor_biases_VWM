function out = mergecell(varargin)
% OUT = MERGECELL(IN1,IN2,...INx,DIM)
%
% MERGECELL simply merges cell arrays (IN1,IN2,...,INx) together by
% appending them along the DIM dimension. There is no requirement that the
% input cell arrays have similar dimensions. The result is a cell array
% large enough to accomodate the sizes of each input array with blanks
% filling the empty spots.
%
% Example:
%
%    a{1} = {'a','b','c','d','e'};
%    a{2} = {'f';'g';'h'};      %-- Notice these are different orientations
%
%    out = mergecell(a{:},1)
%
%     returns:   'a'  'b'  'c'  'd'  'e'  'f'
%                 []   []   []   []   []  'g'
%                 []   []   []   []   []  'h'
%
%
%    out = mergecell(a{:},2)
%
%     returns:   'a'  'b'  'c'  'd'  'e'
%                'f'   []   []   []   []
%                'g'   []   []   []   []
%                'h'   []   []   []   []
%
%
%    out = mergecell(a{1},a{2}',1)
%
%     returns:   'a'  'b'  'c'  'd'  'e'  'f'  'g'  'h'
%
%
% Copyright (C) 2011 Jeremy Brower.
% jeremy.brower@gmail.com
if nargin < 3
    error('MERGECELL needs at least three inputs.');
elseif ~isnumeric(varargin{nargin}) && varargin{nargin} ~= 1 && varargin{nargin} ~=2
    error('The last input argument for MERGECELL needs to be a scalar. Either 1 (merge as rows) or 2 (merge as columns).');
end
for i = 1:nargin-1
    data{1,i} = varargin{i};
end
switch varargin{nargin}
    
    %-- merge into one "row"
    case 1
        totalrows = 0;
        totalcols = 0;
        for i = 1:size(data,2)
            totalrows = max(totalrows,size(data{i},1));
            totalcols = totalcols + size(data{i},2);
        end
        
        out = cell(totalrows,totalcols);
        
        colindex = 1;
        for i = 1:size(data,2)
            for j = 1:size(data{i},2)
                for k = 1:size(data{i},1)
                    out{k,colindex} = data{i}{k,j};
                end
                colindex = colindex + 1;
            end
        end
        
    %-- merge into one "column"
    case 2
        totalrows = 0;
        totalcols = 0;
        for i = 1:size(data,2)
            totalrows = totalrows + size(data{i},1);
            totalcols = max(totalcols,size(data{i},2));
        end
        
        out = cell(totalrows,totalcols);
        
        rowindex = 1;
        for i = 1:size(data,2)
            for j = 1:size(data{i},1)
                for k = 1:size(data{i},2)
                    out{rowindex,k} = data{i}{j,k};
                end
                rowindex = rowindex + 1;
            end
        end
end
