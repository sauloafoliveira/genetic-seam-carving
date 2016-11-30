function [ I, output ] = gcarving( I, params)
    nv = params.nv;
    nh = params.nh;
    
    if ~isfield(params, 'reduceFactor')
        params.reduceFactor = .1;
    end
    
    reduceFactor = params.reduceFactor;
    
    
    if size(I, 3) > 1
        O = rgb2gray(I);
        sp = detectSURFFeatures(O);
        sp = reduce(sp, reduceFactor);
        energyMap = zeros(size(O));
    else
        sp = detectSURFFeatures(I);
        sp = reduce(sp, reduceFactor);
        energyMap = zeros(size(I));
    end
    
    idx = sub2ind(size(O), round(sp.Location(:, 2)), round(sp.Location(:, 1)));
    energyMap(idx) = 1;
    
    output1 = [];
    if nv ~= 0
        itr = abs(nv);
        
        [I, energyMap, output1] = carve(I, energyMap, params, itr);

    end
    
    output2 = [];
    if nh ~= 0
        I = imrotate(I, 90);
        energyMap = imrotate(energyMap, 90);
        
        itr = abs(nh);
        [I, ~, output2] = carve(I, energyMap, params, itr);
        
        I = imrotate(I, -90);

        
    end
    
    
    output = [output1, output2];
end

function [I, energyMap, output] = carve(I, energyMap, params, numberOfSeams)
    count = 0;
    
    feasibleQtd = [];
    Y = [];
    L = [];
    U = [];
    
    generations = [];
    found = [];
    i = 1;
    
    while count < numberOfSeams
        [seams, finalScores, scores, generations(i)] = findseamsGA(energyMap, params);
     
        feasibleQtd(i) = length(finalScores);
        
        smean = mean(scores);
        Y(i) = smean;
        L(i) = smean - min(scores);
        U(i) = max(scores) - smean;
       
        if isempty(seams) 
            continue;
        end
        
        computed_seams = filterseams(seams, numberOfSeams - count);
        
        found(i) = size(computed_seams, 1);
        
        energyMap = carvevseam(energyMap, computed_seams, false);
        I = carvevseam(I, computed_seams, true);
        count = count + found(i);
        
        i = i + 1;
    end
    
    output = {};
    output.cycles = i - 1;
    output.meanFitness = Y;
    output.lowerFitness = L;
    output.upperFitness = U;
    
    output.feasible = feasibleQtd;
    output.processed = found;
    output.generations = generations;
    
end

function [seams, finalScores, scores, generations] = findseamsGA(I, globaParams)
    numberOfGenerations = globaParams.numberOfGenerations;
    populationSize      = globaParams.populationSize;
    alpha = 0;
    
    
    params = {};
    params.numberOfGenerations = numberOfGenerations;
    params.populationSize = populationSize;
    params.geneBound = size(I, 2);
    params.numberOfChromosomes = 1;
    params.numberOfGenes = size(I, 1) * params.numberOfChromosomes;
    params.I = I;
    params.wsize = globaParams.wsize;
    params.CrossoverFraction = .8;
    params.plot = {};
    
    [~, fval, population, scores, generations] = GeneticAlgorithm(params, I, params.wsize);
    
    
    feasible = scores <= (fval + alpha);
    
    fittestIndividuals = population(feasible, :);
    finalScores = scores(feasible);
    % transform into seams
    seams = [];
    for i = 1:size(fittestIndividuals, 1)
        seams = [seams; toseam(fittestIndividuals(i, :))];
    end
end


