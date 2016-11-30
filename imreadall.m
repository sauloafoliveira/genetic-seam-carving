function [ images, filenames ] = imreadall( folder, pattern, varargin )
    images = {};
    
    start = 0;
    
    if nargin > 2
        start = varargin{1};
    end
    
    
    files = dir(folder);
    filenames = {files.name}';
    fun = @(x) (numel(x) > 0) && ((x(1) == 1) || ~start); % garantir so que comece mesmo
   
    indexes = cellfun(fun, strfind(filenames, pattern));
    
    filenames = {files(indexes).name}';
    
    for i = 1:size(filenames, 1)
        filename = [folder, filenames{i}];
        images{i} = imread(filename);
    end

    images = images';
end

