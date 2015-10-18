function d = fpEaglemanDistance(source, target)
%
% Compute Eagleman distance between color matches for a source and target.
% Source and target can be 2D (num letters x 3) or 3D (num subjects x num
% letters x 3). 
%
% d = fpEaglemanDistance(source, target)
%
% output a distance matrix, n subjects in source x n subjects in target
%
% Example:
%   source = fpSimulateData(1,'magnets');
%   target = fpSimulateData(500, 'random');
% 
%   D = fpEaglemanDistance(source, target);


if length(size(source)) < 3
    source = reshape(source, 1, size(source, 2), size(source,3));    
end

if length(size(target)) < 3
    target = reshape(target, 1, size(target, 2), size(target,3));    
end


nsource = size(source, 1);
ntarget = size(target, 1);
nletter = size(source, 2);

if ~isequal(nletter, size(target, 2))
    error('The number of letters in the source (%d) is not the same as the number of letters in the target (%d)',...
        nletter, size(target, 2));
end

d = zeros(nsource, ntarget);

for s = 1:nsource
    for t = 1:ntarget
        tmp = abs(source(s, :,:) - target(t, :,:));
        d(s, t) = sum(tmp(:)) / nletter;
    end
end