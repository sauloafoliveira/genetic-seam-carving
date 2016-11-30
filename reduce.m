function [ reduced ] = reduce(surfPoints, reductionFactor)
    points = surfPoints.Location;
    k = round(reductionFactor * surfPoints.Count);
    [~, centers] = kmeans(points, k, 'MaxIter', 500);
    nearest = dsearchn(points, centers);
    reduced = surfPoints(nearest);
end