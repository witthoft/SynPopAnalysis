% % hierarchical clustering of face similarity data.

load HCPlotExample.mat
% variables we have are
% dists: which are behavioral ratings of the distance between pairs of faces
% faces: cell array of images
% labels:  for each face image


% first we generate our tree from our distances

% clstdef = 'average';%criteria for clustering
clstdef = 'average';

% get the tree structure
tree = linkage(dists, clstdef); %clustering function

% plot the results
figure('name','dendrogram of distances between faces','Color', [1 1 1]);

%set a threshold for giving different groups of branches a color, do this by
%eyeball
ct = 1.5;

% make tree
[H, T, perm] = dendrogram(tree,0,'colorthreshold',ct,'Labels',labels);
% though I haven't done the bookkeeping in this example, you need perm to
% figure out where your labels ended up.  then you will need to resort your
% images
%     [H,T,OUTPERM] = dendrogram(...) generates a dendrogram and returns a
%     vector of the node labels of the leaves shown in the dendrogram,
%     ordered from left to right on a horizontal dendrogram and bottom to
%     top for a vertical dendrogram. OUTPERM is a permutation of the vector
%     1:P, where P is the number of nodes shown



% Add the faces.  you could put them below the x-axis but I just make the
% plot a little bigger to fit them inside the current axes

limits = get(gca,'YLim');
set(gca,'YLim',[-.5 limits(2)]);
hold on;


% set the image size to scale your pictures
% set imsize to be half the image
imsize = .3;

% now draw you images centered in the space under each branch
% again, I didn't show my bookkeeping but you need the perm output from
% dendrogram to keep track of the order of images
for i = 1:length(faces)
    image([i-imsize;i+imsize],[-.25+imsize; ,-.25-imsize],faces{i});
end

% compute cophenetic cluster coefficient (meant to be a measure of the fit
% of the clustering, by comparing distances in tree with original distances

cr = cophenet(tree,dists);

xlabel(char(['cophenetic correlation coefficient = ' num2str(cr,2)]),'FontSize',24,'FontWeight','Bold');
ylabel('distance','FontSize',24);
title(char(['clustering on sim ratings'], ['clusters defined using ' clstdef]),'FontSize',24);


hold off;



set(gca, 'XTickLabel',labels);


set(gca,'FontSize',16);
set(gca,'FontWeight','bold')

% make the lines thicker
set(H(:),'LineWidth',8);
