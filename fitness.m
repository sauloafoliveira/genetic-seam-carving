function [ fv ] = fitness( varargin )
    
   
    population  = varargin{1};
    energyMap   = varargin{2};
    wsize       = varargin{3};
    wby2        = floor(wsize / 2);
    
    fv = zeros(size(population, 1), 1);
    [x, y] = find(energyMap == 1);
    surf_points = [x, y];
   
    dim = size(energyMap);
    
    for i = 1:size(population)
        [seam, fv(i)]  = penalty(population(i, :), dim);
    
        if(fv(i) > 0)
            fv(i) = fv(i) + 1;
            continue;
        end
        
        s = [ flipud(seam)', (1:length(seam))'];
       
        
        
        fv(i) = lambda(surf_points, s, wby2);
    end
    
end


function [seam, pv] = penalty(chr, dim)
    % conserta o range -2,-1,0 -> -1,0,1
    seam = toseam(chr);
    pv = sum(seam < 1) + sum(seam > dim(2));
end



function [fv] = lambda(matchedPoints, seam, wsize)
    fv = 0;

    if numel(seam) ~= 0
        nearestIndexes = dsearchn(seam, matchedPoints);

        nearestPoints = [nearestIndexes, seam(nearestIndexes, 2)];
        N = sqrt(sum(abs(nearestPoints - matchedPoints) .^ 2, 2));


        pointsInside = N < wsize;
        distanceWeight = sum(abs(wsize - N(pointsInside)) ./ wsize);
        fv =  distanceWeight / size(N, 1);
    end
end