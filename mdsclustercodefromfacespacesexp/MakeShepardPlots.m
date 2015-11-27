% make shepard plot
%this script can be run after the mds solutions are arrived at using
%SsRatedAllPairs

%goal is to plot the original distances (similarities as given by subjects
%for example) against the distances in the solution
%good solutions preserve the rank order of the similarities in a linear
%fashion

%need to unpack original similarities so that have a single row
%start by getting uppertriangular portion

%can get a vector of original nonzero distances by finding the indices of
%nonzero elements of the original distance matrix.  use the transpose since
%find seems to search columnwise





% origdists = nonzeros(triu(dists)');
% this doesn't work because identical things actually have a distance of 0,
% so doesn't pull out the right vector
%so need to unpack the upper triangular part of dists

origdists =[]

for i=1:length(dists)
    origdists = [origdists dists(i, i+1:length(dists))]
end




% plots for solutions of all dimensions up to 6
maxdims = 6

n=n+1; figure(n);


whitebg([0 0 0])

for i=1:maxdims
    soldists = pdist(dim_solutions{i});
    subplot(2,3,i);
%     plot(origdists,soldists,'o');
 plot(origdists,soldists,'o');
    hold on;
    xlabel('dissimilarities','FontSize',14)
    ylabel('distances','FontSize',14)
    title([ num2str(i) '-d solution'],'FontSize',14);
    plot(0:.5:6,0:.5:6,'--','color',[1 0 0]);
    axis equal;
end

% test

