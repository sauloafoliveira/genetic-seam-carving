function [x,fval, population,scores, generations] = GeneticAlgorithm(params, I, wsize)
    
    CreationFcn = {@MyCreationTerFn, params.geneBound, params.numberOfChromosomes};
    MutationFcn = {@MyMutationTerFunc, params.geneBound, params.numberOfChromosomes};
    CrossOverFunc = {@MyCrossOverTerFunc, params.numberOfChromosomes};
    FitnessFcn = {@fitness, I, wsize};
    
    options = gaoptimset('PopulationType', 'bitstring', ...
        'Generations',      params.numberOfGenerations, ...
        'PopulationSize',   params.populationSize, ...
        'SelectionFcn',     @selectionroulette, ...
        'CrossoverFcn',     CrossOverFunc, ...        
        'CreationFcn',      CreationFcn, ...
        'MutationFcn',      MutationFcn, ...
        'CrossoverFraction', params.CrossoverFraction, ...
        'Vectorized',       'on', ...
        'PlotFcns',         params.plot);
   
    
    [x, fval, ~, output, population, scores] = ...
                        ga(FitnessFcn, params.numberOfGenes, options);
                    
    generations = output.generations;
end

function [population] = MyCreationTerFn(varargin)
    GenomeLength = varargin{1};
    options = varargin{3};
    geneBound = [1, varargin{4}];
    numberOfChromosomes = varargin{5};

    popSize =  options.PopulationSize;
    population = randi([-2, 0], popSize, GenomeLength);
    
    singleSeamLength = GenomeLength / numberOfChromosomes;
    mask = (0:numberOfChromosomes- 1) * singleSeamLength;
    mask = repmat(mask, popSize, 1);
    tmp  = randi([1, singleSeamLength], popSize, numberOfChromosomes);
    pivotPositions  = tmp + mask;
    
    % set the indexes
    mask = repmat((1:popSize)', 1, numberOfChromosomes);
    idx = sub2ind(size(population), mask, pivotPositions);
    population(idx) = randi(geneBound, popSize, numberOfChromosomes);
   
    
end


%%%%
function [mutatedChildren] = MyMutationTerFunc(parents,~,GenomeLength,...
    ~,~,~,thisPopulation,varargin)
    
    geneBound = varargin{1};
    numberOfChromosomes = varargin{2};
    mutationRate = 0.5;
    mutatedChildren = zeros(length(parents), GenomeLength);
    
    
    individuals = thisPopulation(parents, :);
    for i=1:length(parents)
        chd = internalMutation(individuals(i, :), numberOfChromosomes, mutationRate, geneBound);
        mutatedChildren(i, :) = reshape(chd', 1, []);
    end
end

function [seams] = internalMutation(individual, numberOfChromosomes, mutationRate, geneBound)
     seams = reshape(individual, [], numberOfChromosomes)';
     seamLength = size(seams, 2);
     
     for i = 1:numberOfChromosomes
        
        mutationPoints = find(rand(1, seamLength) < mutationRate);
        [~, pivotPosition] = max(seams(i, :));

        mutationPoints(mutationPoints == pivotPosition) = [];

        mutated = randi([-2, 0], 1, length(mutationPoints));
        seams(i, mutationPoints) = mutated;
        newPivotValue = seams(i, pivotPosition);

        while newPivotValue == seams(i, pivotPosition)
            newPivotValue = randi([1, geneBound]);
        end
        
        seams(i, pivotPosition) = newPivotValue;
     end
end

function [xoverKids] = MyCrossOverTerFunc(parents,~,GenomeLength, ...
    ~,~,thisPopulation, varargin)

    numberOfChromosomes = varargin{1};

    nKids = length(parents)/2;
    xoverKids = zeros(nKids,GenomeLength);
    index = 1;
    
    individuals = thisPopulation(parents, :);
    for i=1:2:nKids
        
        parent1 = individuals(index, :);
        index = index + 1;
        parent2 = individuals(index, :);
        index = index + 1;
        
        [chd1, chd2] = xoverInternal(parent1, parent2, numberOfChromosomes);
        
        xoverKids(i, :)     = chd1;
        xoverKids(i + 1, :) = chd2;
    end
    
    if size(xoverKids, 1) > nKids
        xoverKids = xoverKids(1:nKids, :);
    end
end

function [chd1, chd2] = xoverInternal(parent1, parent2, numberOfChromosomes)

    seams1 = reshape(parent1, numberOfChromosomes, []);
    seams2 = reshape(parent2, numberOfChromosomes, []);
    
    chd1 = [];
    chd2 = [];
    for i=1:numberOfChromosomes
        
        %save pivots
        [v1, pv1] = max(seams1(i, :));
        [v2, pv2] = max(seams2(i, :));
        
        %remove pivots
        s1 = [seams1(1, 1:pv1-1), seams1(1, pv1 + 1:end)];
        s2 = [seams2(1, 1:pv2-1), seams2(1, pv2 + 1:end)];
        
        % choose two (nonequal) crossover points
        sz  = length(s1) - 1;
        xOverPoint1 = ceil(sz * rand);
        xOverPoint2 = ceil(sz * rand);
        while(xOverPoint2 == xOverPoint1)
            xOverPoint2 = ceil(sz * rand);
        end

        % Deal with the case where the splice wraps around the ends.
        if(xOverPoint1 < xOverPoint2)
            left = xOverPoint1;
            right = xOverPoint2;
        else
            left = xOverPoint2;
            right = xOverPoint1;
            swap = parent1;
            parent1 = parent2;
            parent2 = swap;
        end
        
        
        string1 = [ s1(1:left), s2((left + 1):right), s1((right + 1):end)];
        string2 = [ s2(1:left), s1((left + 1):right), s2((right + 1):end)];
        
        %insert pivot
        
        chd1 = [chd1, string1(1:pv1-1), v1, string1(pv1:end)];
        chd2 = [chd2, string2(1:pv2-1), v2, string2(pv2:end)];
        
        
        
    end
end