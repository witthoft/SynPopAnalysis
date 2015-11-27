% hierarchical clustering of face similarity data.
% takes nx3 matrix of data as input where columns are face 1, face 2,
% similarity rating.  may set flags for computing on subsets of data and
% different definitions for determining a cluster

% set path, variables, and load data


clear
whitebg([1 1 1]);

% figure counter
n = 1;

% Number of subjects.
subjects = 1;

% number of dimensions in the solution
soldims = 2;


% Put subject matrix names into a cell structure.
% create cell array
ratings = cell(subjects,1);

thepath = '/Users/nathan/Experiments/FaceSpaces/ResemblancesWorking/resviewsresults';

% this is a hack to get the angles info. when you pull out the ids they get
% put into numerical order so can assume that this vector can be appended
% to array ids
% specify angle of each stimulus
angles = [0 30 60 90 0 30 60 90 0 30 60 90 0 30 60 90 0 30 60 90 0 30 60 90];
compare = [0 30 60 90];

% 2.08 data
% cd /Users/nathan/Experiments/FaceSpaces/ResemblancesWorking/resviewsresults/2.08.08results/mds
% load 'twoEightSubs.mat'
% ids = unique(Alisha(:,1));
% load '2.08.pixsims.mat';
% ratings{1}=PixSimData; subname='Pixel Differences';


% 2.25 data
% cd /Users/nathan/Experiments/FaceSpaces/ResemblancesWorking/resviewsresults/2.25results/mds
% load '2.25data.mat';
% ids = unique(AllSubs(:,1));
% % %
% %  need to fix, this data is still incorrect
% % ratings{1}=AllSubs; subname = 'allSubs';
% % ratings{1}=Amanda; subname = 'Amanda';
% % ratings{1}=Andrew; subname = 'Andrew';
% % ratings{1}=Aryana; subname = 'Aryana';
% % ratings{1}=David; subname = 'David';
% % ratings{1}=David2; subname = 'David2';
% % ratings{1}=Jeremy;
% % ratings{1}=NoJeremorD2; subname = 'Average wo J and D2';
% load '2.25PixSimData';
% ratings{1}=PixSimData;  subname = 'Pixel Differences';

% 3.24 tri data
cd ~/Experiments/FaceSpaces/ResemblancesWorking/resviewsresults/3.24.results/tri/
load '3.24ratings.mat'
ids = unique(allsubst(:,1));
% ratings{1}=chloet; subname = 'Chloe';
% ratings{1}=dariust; subname = 'Darius';
% % ratings{1}=lestert; subname = 'Lester';
ratings{1}=shinet; subname ='Shine';
% ratings{1}=yongt; subname = 'Yong';
% ratings{1}=allsubst; subname = 'Average';
% % ratings{1}=allnoshine; subname = 'Average wo Shine';
% load '3.24pixsims.mat';
% ratings{1}=PixSimData; subname = 'Pixel Differences';

ct = 2.2 %color threshhold

% 3.11.data
% cd /Users/nathan/Experiments/FaceSpaces/ResemblancesWorking/resviewsresults/3.11.results/mds;
% load '311SubData.mat';
% ids = unique(AllSs(:,1));
% ratings{1}=Adrian; subname = 'Adrian';
% ratings{1}=AllSs;  subname = 'Average';
% ratings{1}=Ideal; subname = 'Ideal';
% ratings{1}=Joel; subname = 'Joel';
% ratings{1}=Darius; subname = 'Darius';
% ratings{1}=Beth; subname = 'Beth';
% ratings{1}=AllNoBeth; subname = 'Average wo Beth';
% ratings{1}=Ideal;

% 4.15 tri data

% cd /Users/nathan/Experiments/FaceSpaces/ResemblancesWorking/resviewsresults/4.15results
% load '415data.mat';
% ids = unique(allsubs(:,1))
% % ratings{1}=allsubs; subname = 'average';
% % ratings{1}=erikka; subname = 'erikka';
% % ratings{1}=greg; subname = 'greg';
% % ratings{1}=steph; subname = 'steph';
% load '4.15pixsims.mat';
% ratings{1}=PixSimData; subname = 'Pixel Differences';


% add column specifiying angles
ids = [ids [1:length(ids)]' angles'];

% pull out just the desired angles
% first find right ids and angles
subset = [];
for i = 1:length(compare)
    subset = [subset; ids( find(ids(:,3)==compare(i)),:,:)]; %looks for angles that match those specified in compare and gets just those
end
ids = [sort(subset) [1:length(subset)]']; %sort them in ascending order
% then take out just the desired subset of data
% a little tricky.  can start by finding the values in the first column of
% ratings that correspond to the faces you want
subset=[];

for i=1:length(ids)
    subset=[subset; ratings{1}(find(ratings{1}(:,1)==ids(i)),:)]
end

% then get the ones that have the right value in the second column
subset2=[];
for i=1:length(ids)
    subset2=[subset2; subset(find(subset(:,2)==ids(i)),:)]
end

% need to remove comparisons of stim to itself
% find places where two faces are not the same

dratings=subset2(find(subset2(:,1)~=subset2(:,2)),:)

% now need to reshuffle the order of the similarities to be in a way that
% linkage or cluster expects.  both take the output of pdist which treats
% each row of a matrix as an item and each column as an observation (or
% distance along a dimension).  the output of pdist is a single row vector
% which has the distances between pairs of objects.  these are ordered such
% that the first item is compared to the rest, then the second, then the
% third .... so if there are 4 items in a, then pdist(a) will have (n-1)!
% entries, with the comparsions being [1 to 2, 1 to 3, 1 to 4, 2 to 3,...3
% to 4].  our data is already in distances, but needs to be resorted so
% that it has the same form as if it had passed through pdist instead of
% people


% unfortunately, the right number might appear in either column of dratings
% so need to go through and pick out the elements. ids is sorted, so you
% can go through it in order to get things into the right columns.

oratings=[];
oaltnames = cell(length(ids),1);

for i=1:length(ids)
    %     find ids in column 1 equal to the ith values of ids
    index = find(dratings(:,1)==ids(i));
    oratings = [oratings; dratings(index,:)];
    %     then want to remove these rows so they won't be the reselected
    dratings(index,:)=[];
    %     now get the ids in column equal to the ith value of ids
    index = find(dratings(:,2)==ids(i));
    %      need to swap first and second column
    oratings = [oratings; dratings(index,2) dratings(index,1) dratings(index,3)];
    %     remove these rows
    dratings(index,:)=[];


end


% get our similarities

sims = oratings(:,3)';

% convert to distances

dists = abs(sims-7);
figure(n); hist(dists); n=n+1;

% now do clustering

clstdef = 'average';%criteria for clustering

tree = linkage(dists, clstdef); %clustering function

% plot the results

figure(n); n=n+1;
% if you want the labels to be the same as the image names
[H, T, perm] = dendrogram(tree,0,'colorthreshold',ct,'Labels',num2str(ids(:,1)));

labels = get(gca,'XTickLabels')



% add faces.
% if I look under get(gca,'YTickLabel') that returns the list of label
% names which is fortunately the same as the name of the face. moreover
% since their position on the y-axis is the same as their index in the
% array its easy to determine the location in y-coordinates.
% two possibilities for the x coordinates.  one is to set the y axis to
% cross at some negative value or eliminate it and just extend the x axis
% to negative values.  both have the effect of opening up whitespace to
% write the images to.

limits = get(gca,'YLim');
set(gca,'YLim',[-.5 limits(2)]);
hold on;
% made space now just need to write each image to be .5 on each side, with
% the left boundary at 0, the right boundery at -.5, the top at i+.25 the
% bottom at i-.25 where i:1 num stims
%  load faces  need those in order on axis so



faces = cell(1,length(labels));
for ctr = 1:length(labels)
    pic = char([strcat(strtrim(labels(ctr)),'.jpg')])
    faces{ctr} = imread(pic);
end


% H = dendrogram(tree,0,'colorthreshold',2.1,'Labels',altnames);

% set imsize to be half the image
imsize = .3;

% now add images
image([1-imsize;1+imsize],[-.25+imsize; ,-.25-imsize],faces{1});

%now the rest  need to paramterize this so that it only uses as many faces
%as in subset
for i = 2:length(ids)
    image([i-imsize;i+imsize],[-.25+imsize; ,-.25-imsize],faces{i});
end
% compute cophenetic cluster coefficient (meant to be a measure of the fit
% of the clustering, by comparing distances in tree with original distances

cr = cophenet(tree,dists);

xlabel(char(['cophenetic correlation coefficient = ' num2str(cr,2)]),'FontSize',24,'FontWeight','Bold');
ylabel('distance','FontSize',24);
title(char(['clustering on ' subname ' ratings'], ['clusters defined using ' clstdef]),'FontSize',24);


hold off;

% use cluster to find clusters.  couple of ways of defining thresholds.  we
% can use arbitrary height from looking at the tree

% T=cluster(tree,'Cutoff',4,'Criterion','Distance');
% for i=1:max(T)
%     disp(['cluster ' num2str(i) ' is']);
%     ind = find(T==i);
%     for j=1:length(ind)
%         disp(ids(ind(j))) %the indices spat out by cluster refer back to the original items through the output of linkage
%     end
% end



% if you want the labels to be easier to read and show the structure of the
% morphs  fortunately, when we pulled out the ids we wanted to use we saved
% their original indices in the second column.  since ids and altnames have
% the same order, we just get the altnames that have the same indices as
% the ids

% first check that the altnames exist
if exist('altnames')
    new_labels = cell(length(ids),1)
    for i=1:length(ids)
        new_labels{i}=altnames{ids(i,2)};
    end


    % now have the right set of names but need to resort them to be in the same
    % order as the labels.  the order for the labels is given to perm by the
    % function dendrogram

    for i=1:length(new_labels)
        labels{i}=new_labels{perm(i)}
    end


    set(gca, 'XTickLabels',labels);
end

set(gca,'FontSize',16);
set(gca,'FontWeight','bold')

% make the lines thicker
set(H(:),'LineWidth',8);
