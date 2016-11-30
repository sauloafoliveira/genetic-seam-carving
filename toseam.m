function [ seams ] = toseam( chr, varargin )
    if nargin > 1
        numberOfChromosomes = varargin{1};
    else
        numberOfChromosomes = 1;
    end
    
    individuals = reshape(chr, [], numberOfChromosomes)';
    seams = [];
    for i = 1:numberOfChromosomes
        seams = [seams; internal_toseam( individuals(i, :) )];
    end
end

function [ seam ] = internal_toseam( chr )
    [value, pivot] = max(chr);
    chr = chr + 1;
    chr(pivot) = value;
    
    seam(pivot) = chr(pivot);
    for i = (pivot + 1):length(chr)
        seam(i) = seam(i - 1) + chr(i);
    end
    
    seam(pivot) = chr(pivot);
    
    for i = (pivot-1):-1:1
        seam(i) = seam(i + 1) + chr(i);
    end
    seam(pivot) = chr(pivot);
    seam = round(seam);
end

