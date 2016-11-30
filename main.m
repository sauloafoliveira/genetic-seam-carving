clear all; close all; clc;

I = imread('bicycle1.jpg');
dim = size(I);

params = {};
params.nv = -100;
params.nh = 0;
params.wsize = 27;
params.reduceFactor = .10;
params.numberOfGenerations = 80;
params.populationSize = 40;

[O, output] = geneticseamcarving(I, params);


imwrite(O, 'sGSC_bicycle.jpg');

imshow([I O]);
title(sprintf('%d cycles and %d generations (avg)', ...
                output.cycles, mean(output.generations)), 'FontSize', 14);


%% output var

% * output.cycle => number of ga runs
% * output.[meanFitness, lowerFitness, upperFitness] => fitness measurements
% for each cycle
% * output.feaible => number of valid seams detected at each cycle.
% * output.processed => number of  seams processed at each cycle.
% * output.generations => number of generations at each cycle.
