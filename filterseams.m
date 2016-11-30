function [ uncrossSeams ] = filterseams( seams, count )
   
    [ uncrossSeams ] = internal_filterseams( seams, count );
end

function [ uncrossSeams ] = internal_filterseams( seams, count )
    % remove those that intersect each other
    hasCrossingIndividuals = 1;
    
    cellseams = num2cell(seams, 2);
    while hasCrossingIndividuals
        hasCrossingIndividuals = 0;
        
        collisions = zeros(size(seams, 1), 1);
        sequence = 1:size(seams);
        
        for i = sequence
            y = seams(i, :);
            t  = cellfun(@(x)(sum((x - y) == 0)), ...
                cellseams(sequence ~= i, :), 'UniformOutput', false);
            collisions(i) = sum(cell2mat(t));
        end
        
        if  max(collisions) == 0
            break;
        end
        
        
        if size(seams, 1) == 1
            break;
        end
        
        [~, indexes] = sort(collisions);
        theOneWithMoreColisions = indexes(end);
        
        seams(theOneWithMoreColisions, :) = [];
        cellseams(theOneWithMoreColisions, :) = [];
        hasCrossingIndividuals = 1;
    end
    
    if size(seams, 1) > count
        uncrossSeams = seams(1:count, :);
    else
        uncrossSeams = seams;
    end
end

