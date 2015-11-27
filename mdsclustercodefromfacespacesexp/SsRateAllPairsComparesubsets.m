% Multidimensional scaling for similarity data. Very ugly.  
%For this particular version subjects rated every possible pair of 6 faces
%at 4 views which is 24!/2!(24-2)! = 276 pairs, + every stimulus with
%itself for a total of 300 pairwise distances.   This when finished should
%have complementary plots for both the individual and average data.
%want the cross-correlation matrices

% goal of this version is to allow specification of which subsets of data
% to look at.  this could be done by just inputting those at the data level
% but seems nicer to allow choice rather than making all possible pairs


%since the subsets are views (0, 30, 60, 90) and the names of the images
%end in 1,3,6 or 9 should be easy to pull just those out.  trick will be
%fixing the plotting

% want the mds solution for each subject and average, as well as subjects
% solutions rotated onto one another
% scree plots to see what the best dimensional explanation might be
% 

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

% decide what angles to compare, must be at least two.  for fewer angles
% may want to reduce the dimensionality of solution?

compare = [0 30 60 90];

%fill cells with either average data or individual data

% 2.25 data
% cd /Users/nathan/Experiments/FaceSpaces/ResemblancesWorking/resviewsresults/2.25results/mds
% load '2.25data.mat';
% ids = unique(AllSubs(:,1));
% % 
%  need to fix, this data is still incorrect 
% ratings{1}=AllSubs; subname = 'allSubs';
% ratings{1}=Amanda; subname = 'Amanda';
% ratings{1}=Andrew; subname = 'Andrew';
% ratings{1}=Aryana; subname = 'Aryana';
% ratings{1}=David; subname = 'David';
% ratings{1}=David2; subname = 'David2';
% ratings{1}=Jeremy;
% ratings{1}=NoJeremorD2; subname = 'Average wo J and D2';

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

% 2.08 data

% cd /Users/nathan/Experiments/FaceSpaces/ResemblancesWorking/resviewsresults/2.08.08results/mds
% load 'twoEightSubs.mat'
% ids = unique(Alisha(:,1));

% ratings{1}=Alisha; subname = 'Alisha';
% ratings{1}=Chee; subname = 'Chee';
% ratings{1}=Salina; subname = 'Salina';
% ratings{1}=Sheila; subname = 'Sheila';
% ratings{1}=Jason; subname = 'Jason';
% ratings{1}=du; subname = 'Du';
% ratings{1}=cora; subname = 'Cora';
% ratings{1}=Alice; subname = 'Alice';
% ratings{1}=SSRatedAllPairs; subname = 'Average';
% ratings{1} = AllSubs; subname = 'AllSubs';


% 3.24 tri data
% cd ~/Experiments/FaceSpaces/ResemblancesWorking/resviewsresults/3.24.results/tri/
% load '3.24ratings.mat'
% ids = unique(allsubst(:,1));
% ratings{1}=chloet; subname = 'Chloe';
% ratings{1}=dariust; subname = 'Darius';
% ratings{1}=lestert; subname = 'Lester';
% ratings{1}=shinet; subname ='Shine';
% ratings{1}=yongt; subname = 'Yong';
% ratings{1}=allsubst; subname = 'Average';
% ratings{1}=allnoshine; subname = 'Average wo Shine';
% load '3.24pixsims.mat';
% ratings{1}=PixSimData; subname = 'Pixel Differences';
ls


%4.1 tri data
%  this one turned out not so great.  may have errors in the data so need
%  to double check
% cd ~/Experiments/FaceSpaces/ResemblancesWorking/resviewsresults/4.1results/tri/
% load '4.1tridata.mat'; 
% ids = unique(allsubs(:,1));
% ratings{1}=allsubs; subname = 'Average';
% need to load all the subjects data
% ratings{1}=allsubsnoted; subname = 'Average wo Ted';
% ratings{1}=brandon; subname = 'Brandon';
% ratings{1}=garrett; subname = 'Garrett';
% ratings{1}=holly; subname = 'Holly';
% ratings{1}=su; subname = 'Su';
% ratings{1}=valeri; subname = 'Valeri';
% ratings{1}=teddy; subname = 'Teddy';


% 4.8 rect data

% cd /Users/nathan/Experiments/FaceSpaces/ResemblancesWorking/resviewsresults/4.8results
% load '4.8data.mat';
% ids = unique(allsubs(:,1));
% % 
% ratings{1}=allsubs; subname = 'average';
% ratings{1}=jonathan; subname = 'jonathan';
% ratings{1}=jin; subname = 'jin';
% ratings{1}=maria; subname='maria';
% ratings{1}=carlo; subname='carlo';
% ratings{1}=charlie; subname='charlie';
% ratings{1}=nina; subname='nina';
% load '4.28pixsimsdata.mat';
% ratings{1}=PixSimData;

% 5.15 tri data

% cd /Users/nathan/Experiments/FaceSpaces/ResemblancesWorking/resviewsresults/4.15results
% load '415data.mat';
% ids = unique(allsubs(:,1))
% ratings{1}=allsubs; subname = 'average';
% ratings{1}=erikka; subname = 'erikka';
% ratings{1}=greg; subname = 'greg';
% ratings{1}=steph; subname = 'steph';
% load '4.15pixsims.mat';
% ratings{1}=PixSimData; subname = 'Pixel Differences';

% 8.18 tri data
cd /Users/nathan/Experiments/FaceSpaces/ResemblancesWorking/resviewsresults/8.18results
load '8.18data.mat';
% ids = unique(allsubs(:,1))
% ratings{1}=allsubs; subname ='average';
% ratings{1}=allsubsnoCorN; subname ='average no Chris or Nichole';
ids = unique(allsubsnoCNMa(:,1))
ratings{1}=allsubsnoCNMa; subname ='average no CNMa';
% ratings{1}=chris; subname = 'chris'  % this person behaved randomly
% ratings{1}=justin; subname = 'justin'  %ok 
% ratings{1}=katerina; subname = 'katerina'  %very good
% ratings{1}=marina; subname = 'marina' %pretty bad as expected
% ratings{1}=matt; subname = 'matt'  %very good
% ratings{1}=michelle; subname = 'michelle' %ok  mixed up 90s
% ratings{1}=nichole; subname = 'nichole' %pretty random
% 
% load 'PixSims.mat'
% ratings{1}=PixSimData; subname = 'Pixel Differences'


% 8.18far  tri data
% cd /Users/nathan/Experiments/FaceSpaces/ResemblancesWorking/resviewsresults/8.18results
% load '8.18fardata.mat';
% ids = unique(allsubs(:,1))
% % ratings{1}=allsubs; subname ='average';
% % ratings{1}=bill; subname = 'bill';
% % ratings{1}=donald; subname = 'donald';
% ratings{1}=junkyu; subname = 'junkyu';


% Get id values for the face numbers.  

%add second column which is counter (i.e. 1 to number of unique faces)
% add third column which specifies viewpoint angle
ids = [ids [1:length(ids)]' angles']

% Construct similarity matrices.
% Make cell array for similarities of same size as number of subjects
similarities = cell(subjects,1);
%Make cell array for distances
distances = cell(subjects,1);




% replace 0s (no rating) with NANs

ratings{1}(find(ratings{1}(:,3)==0),3)=NaN;

% plot similrities and distances


all_sims = [];


% make similarity matrices, 1 for each subject
% need to ammend this to only choose subset asked for.
% one way is to just reduce ids to subset of angles and then go from there
% leaving reast of code untouched (excellent)!
subset = [];

for i = 1:length(compare)
subset = [subset; ids( find(ids(:,3)==compare(i)),:,:)]
end

ids = [sort(subset) [1:length(subset)]']


% may want to have alternate labels (besides the image names) which will be
% put in a variable altnames loaded with the data.  it is a cell array with
% the same number of names as images, and sorted to correspond with ids.
% need to just select the subset of these we want, and ids(:,2) has the
% original indices
if exist('altnames')
    new_labels = cell(length(ids),1);
    for i=1:length(ids)
        new_labels{i}=altnames{ids(i,2)}
    end
end


for subject = 1:subjects
    %makes square matrix of size n x n where n is the number of unique faces
    sims = zeros(length(ids),length(ids));

    %get data for 1 subject each data file has 3 columns, face1, face2, and
    %the similarity rating
    data = ratings{subject};
    %number of rows and columns
    [di dj] = size(data);
    
    
    %go to each row
    for i = 1:di
        %fill sims array (like an xcel pivot table)
        %basically for a given row in the data file, want to find the corresponding point in the nxn matrix where the rating goes
        %ids has two columns,ids(:,1) has original face database number, and ids(:,2) has new numbering from 1 to the number of unique faces
        %data(:,1) has face1, data(:,2) has face2, data(:,3) has rating
        %go to row i in data and         
        %find row of ids in which face1 is equal to shown face1
        %then find row of ids in which face is equal to shown face2
        %this gives you sims(face1, face2) which is set equal to the rating]
        
        


        
        sims(find(ids(:,1)==data(i,1)),find(ids(:,1)==data(i,2))) = data(i,3);
    end

    %        this only works if the ordering is correct in the data  need it to
    %        sort through all the data and reorder it in a sensible way
    %         ids(1,:) has the faces  maybe easiest to resort data  could
    %         just run it in the other direction....
        
     for i = 1:di
        %fill sims array (like an xcel pivot table)
        %basically for a given row in the data file, want to find the corresponding point in the nxn matrix where the rating goes
        %ids has two columns,ids(:,1) has original face database number, and ids(:,2) has new numbering from 1 to the number of unique faces
        %data(:,1) has face1, data(:,2) has face2, data(:,3) has rating
        %go to row i in data and         
        %find row of ids in which face1 is equal to shown face1
        %then find row of ids in which face is equal to shown face2
        %this gives you sims(face1, face2) which is set equal to the rating]
        
        


        
        sims(find(ids(:,1)==data(i,2)),find(ids(:,1)==data(i,1))) = data(i,3);
    end
    
    
    
    % Draw histogram of similarity ratings.
    %  want to add figure headings
    figure(n);  n=n+1;
    hist(sims(:),7);
    
    xlabel('rating','FontSize',24);ylabel('number of pairs','FontSize',24);title([subname ' ratings'],'FontSize',24);
    
    
    %start with an empty matrix, and on each pass add another colum
    %corresponding  to a single subject's ratings
%     sims=fliplr(flipud(sims));
    all_sims = [all_sims sims(:)];
    % then Turn similarities into distances.
    dists = max(max(sims))- sims;  %rating is from 1 to 7
    %fill whole matrix with upper triangular part (makes sure its symmetric)
    dists = triu(dists) + triu(dists)';
    distances{subject} = dists;
    distances{subject};
    %should end up with two arrays, sims and a cell array dists
end

% Output correlation matrix for subjects' ratings.
%i.e. how correlated subjects are with one another
%this has a problem in that when subjects all do not do some comparisons
%(say across views) then there will be 0s in the matrix
%one possibility is to look for them and replace them with nans
adj_sims = all_sims;




% should shift this up to show each subjects matrix
% figure(n); 
% subplot(3,1,1); image(triu(sims)); colormap(bone(8)); xlabel('similarities','Color',[1,0,0],'FontSize',18); colorbar;
% subplot(3,1,2); image(triu(dists)); colormap(bone(8)); xlabel('distances','Color',[0,0,0],'FontSize',18); colorbar;



% make empty square matrix
sum_distance_matrix = zeros(length(ids));


%if more than one subject then    
% Calculate mean distance matrix.

if subjects>1
    %sum each subject's data
    for subject = 1:subjects
        sum_distance_matrix = sum_distance_matrix + distances{subject};
    end
    %divide by number of subjects
    sum_distance_matrix;
    subjects;
    mean_distance_matrix = sum_distance_matrix ./ subjects
else
    mean_distance_matrix = distances{subject};
end
%cheat to set diagonal values to zero (i.e. distance between identical
%items)
for i=1:size(max(ids))
    mean_distance_matrix(i,i)=0;
end

% %plot mean distance matrix
% subplot(3,1,3); image(triu(mean_distance_matrix)); colormap(bone(8)); xlabel('distances','Color',[0,0,0],'FontSize',18); colorbar;
%  


% large plot of similarities
n=n+1;
figure(n);
imagesc(triu(sims)); colormap(bone(8)); xlabel('similarities','Color',[1,0,0],'FontSize',18); colorbar;

n=n+1;

% Get scaling solution for mean distance matrix.  this just gets one for a
% chosen dimensionality 
opts = statset('MaxIter',2000)
[mean_mds_solution,stress,disparities] = mdscale(mean_distance_matrix,soldims,'Criterion','stress','Start','random','Options',opts);






% Load face images.
%gives you a cell array with each element corresponding to a picture of a face

faces = cell(1,length(ids));
for ctr = 1:length(ids)
    file = [num2str(ids(ctr,1)) '.jpg']
    faces{ctr} = imread(file);
end


% Some plotting parameters.
imsize = 0.3;
fig_posn = [0 0 1000 1000];


% do all possible dimensionalities and plot stress measures (scree plot) to
% see underlying dimensionality of data

%highest dimensionality is number of faces
numcomparisons = length(ids)/4;

%make matrix for stresses  first column is stresses, second is counter
strs = [zeros(numcomparisons,1), [1:numcomparisons]'];

%make cell array for mean solutions 
dim_solutions = cell(numcomparisons);
%make cell array for disparities
dim_disparities = cell(numcomparisons);

%then iterate through solutions and get stresses

for i =1:numcomparisons
    %do mds with dimensionality i
    opts = statset('MaxIter',2000)
    %would be better to write these to a series of cell arrays
    [mds_solution,stress,disparities] = mdscale(mean_distance_matrix,i,'Criterion','stress','Start','random','Options',opts);
    %compute
    strs(i,1)=stress;
    dim_solutions{i}= mds_solution;
    dim_disparities{i}= disparities;
end






% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
% plot stress x dimensionality
figure(n);
n=n+1;
plot(strs(:,2),strs(:,1));
xlabel('number of dimensions in solution','FontSize',18);
ylabel('stress','FontSize',18);
title('goodness of fit for solutions of different dimensionality');








% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 


% want to plot arbitrary pair of dimensions in a figure
% need to know dimensionality of the solution
% and which two dimensions to plot

num_dims_insol = 2;
dim1 = 1;
dim2 = 2;


%get solution of right dimensionality
sol_to_plot = dim_solutions{num_dims_insol};

% fill y dim of plot if solution is only 1 dimensional
if num_dims_insol == 1
    sol_to_plot(:,2)=0;
end


%now plot chosen dimensions



imsize = 0.3;
fig_posn = [0 0 1000 1000];

n=n+1;
F=figure(n); 


%plot first face
image([sol_to_plot(1,dim1)-imsize;sol_to_plot(1,dim1)+imsize],[sol_to_plot(1,dim2)-imsize; sol_to_plot(1,dim2)+imsize],faces{1});
hold on;

%now the rest  need to paramterize this so that it only uses as many faces
%as in subset
for ctr = 2:length(ids)
   image([sol_to_plot(ctr,dim1)-imsize; sol_to_plot(ctr,dim1)+imsize],[sol_to_plot(ctr,dim2)-imsize; sol_to_plot(ctr,dim2)+imsize],faces{ctr}); 
end


% will always be 6 identities
nv = length(compare);

% colors for plot
labelmap = [1 0 0; .5 0 .5; 0 0 1; 1 .5 0; 1 1 0; 0 1 0];


% for i = 1:length(ids)/nv
%     text(sol_to_plot(i*nv-(nv-1):i*nv,dim1),sol_to_plot(i*nv-(nv-1):i*nv,dim2),num2str([i*nv-(nv-1):i*nv]'),'Color',labelmap(i,:), 'FontSize',18);
% end



% for i = 1:length(ids)/nv
%     text(sol_to_plot(i*nv-(nv-1):i*nv,dim1),sol_to_plot(i*nv-(nv-1):i*nv,dim2) ,num2str(ids(i*nv-(nv-1):i*nv)'),'Color',labelmap(i,:), 'FontSize',24,'FontWeight','Bold');
% end


% alternatively plot names from morph space (supplied by variable altname
% that would be loaded with the data
for i = 1:length(ids)/nv
    text(sol_to_plot(i*nv-(nv-1):i*nv,dim1),sol_to_plot(i*nv-(nv-1):i*nv,dim2) ,new_labels(i*nv-(nv-1):i*nv)','Color',labelmap(i,:), 'FontSize',20,'FontWeight','Bold');
end

xlabel(['dimension' num2str(dim1)],'Color',[1,0,0],'FontSize',18);  ylabel(['dimension' num2str(dim2)],'Color',[0,1,0],'FontSize',18);
title(['mds on ' subname ' ratings : '    num2str(num_dims_insol) ' dimensional solution'],'FontSize',24);

%find how big the plot needs to be
min_x = min(sol_to_plot(:,dim1));
max_x = max(sol_to_plot(:,dim1));
min_y = min(sol_to_plot(:,dim2));
max_y = max(sol_to_plot(:,dim2));

vals = abs([min_x max_x min_y max_y]);

%largest distance + some wiggle room
val = max(vals) + imsize;
%set axes
axis([-val val -val val]);
colormap([[0:1/255:1]', [0:1/255:1]',[0:1/255:1]'])
set(F,'Position',fig_posn);
hold off;





% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
% WANT TO PLOT ARBITRARY PAIR OF DIMENSIONS IN A 3DFIGURE
% NEED TO KNOW DIMENSIONALITY OF THE SOLUTION
% AND WHICH TWO DIMENSIONS TO PLOT

num_dims_insol = 3;
dim1 = 1;
dim2 = 2;
dim3 = 3;
%get solution of right dimensionality
sol_to_plot = dim_solutions{num_dims_insol};

%now plot chosen dimensions



imsize = 0.3;
fig_posn = [0 0 1000 1000];

n=n+1;
F=figure(n); 


%plot first face
image([sol_to_plot(1,dim1)-imsize;sol_to_plot(1,dim1)+imsize;],[sol_to_plot(1,dim2)-imsize; sol_to_plot(1,dim2)+imsize],faces{1});
hold on;

%now the rest
for ctr = 2:length(ids)
   image([sol_to_plot(ctr,dim1)-imsize; sol_to_plot(ctr,dim1)+imsize],[sol_to_plot(ctr,dim2)-imsize; sol_to_plot(ctr,dim2)+imsize],faces{ctr}); 
end



% plot names of images on screen
% for i = 1:length(ids)/nv
%     text(sol_to_plot(i*nv-(nv-1):i*nv,dim1),sol_to_plot(i*nv-(nv-1):i*nv,dim2),sol_to_plot(i*nv-(nv-1):i*nv,dim3) ,num2str(ids(i*nv-(nv-1):i*nv)'),'Color',labelmap(i,:), 'FontSize',18);
% end

% alternatively plot names from morph space (supplied by variable altname
% that would be loaded with the data
for i = 1:length(ids)/nv
    text(sol_to_plot(i*nv-(nv-1):i*nv,dim1),sol_to_plot(i*nv-(nv-1):i*nv,dim2),sol_to_plot(i*nv-(nv-1):i*nv,dim3) ,new_labels(i*nv-(nv-1):i*nv)','Color',labelmap(i,:), 'FontSize',20,'FontWeight','Bold');
end

% label axes
xlabel(['dimension' num2str(dim1)],'Color',[1,0,0],'FontSize',18);  ylabel(['dimension' num2str(dim2)],'Color',[0,0,1],'FontSize',18);
 zlabel(['dimension' num2str(dim3)],'Color',[0,1,0],'FontSize',18);
 
title(['mds on ' subname ' ratings : '   num2str(num_dims_insol) ' dimensional solution'],'FontSize',24);

%find how big the plot needs to be
min_x = min(sol_to_plot(:,dim1));
max_x = max(sol_to_plot(:,dim1));
min_y = min(sol_to_plot(:,dim2));
max_y = max(sol_to_plot(:,dim2));
min_z = min(sol_to_plot(:,dim2));
max_z = min(sol_to_plot(:,dim3));

vals = abs([min_x max_x min_y max_y min_z max_z]);

%largest distance + some wiggle room
val = max(vals) + imsize+1;
%set axes
axis([-val val -val val -val val]);
colormap([[0:1/255:1]', [0:1/255:1]',[0:1/255:1]'])
set(F,'Position',fig_posn);
hold off;





cd (thepath);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 

%plot points surrounded by convex hull.  should allow visualization of
%groupings more easily

% first need num dimensions in solution and dimensions to plot
n=n+1;
figure(n);


num_dims_insol = 3;
dim1 = 1;
dim2 = 2;
dim3 = 3;
%get solution of right dimensionality
sol_to_plot = dim_solutions{num_dims_insol};

% set up colormap for figure to get colors right

colormap([1 0 0; .5 0 .5; 0 0 1; 1 .5 0; 1 1 0; 0 1 0]);

    caxis([1,6]);   

col = [1 1 1 1
        2 2 2 2
        3 3 3 3 
        4 4 4 4
        5 5 5 5
        6 6 6 6];

%now find convex hull for each of the identities

% 
% hull=[sol_to_plot(:,1),sol_to_plot(:,2),sol_to_plot(:,3)];
% for i = 1:length(ids)/nv
%      patch(sol_to_plot(i*nv-(nv-1):i*nv,1),sol_to_plot(i*nv-(nv-1):i*nv,2),sol_to_plot(i*nv-(nv-1):i*nv,3),labelmap(i,:));alpha(.5);
% %     trisurf(hull,sol_to_plot(i*nv-(nv-1):i*nv,1),sol_to_plot(i*nv-(nv-1):i*nv,2),sol_to_plot(i*nv-(nv-1):i*nv,3),col(i,:)); alpha(.5);
% end

hull1 = convhulln(sol_to_plot(1:4,:));
trisurf(hull1,sol_to_plot(1:4,1),sol_to_plot(1:4,2),sol_to_plot(1:4,3),col(1,:)); alpha(.5);
% fill3(hull1(:,1),hull1(:,2),hull1(:,3),col1');
hold on;
% patch(hull1(2,1),hull1(2,2),hull1(2,3),[1 0 0]);
% patch(hull1(3,1),hull1(3,2),hull1(3,3),[1 0 0]);
% patch(hull1(4,1),hull1(4,2),hull1(4,3),[1 0 0]);



hull2 = convhulln(sol_to_plot(5:8,:));
trisurf(hull2,sol_to_plot(5:8,1),sol_to_plot(5:8,2),sol_to_plot(5:8,3),col(2,:));alpha(.5);
hull3 = convhulln(sol_to_plot(9:12,:));
trisurf(hull3,sol_to_plot(9:12,1),sol_to_plot(9:12,2),sol_to_plot(9:12,3),col(3,:)); alpha(.5);
hull4 = convhulln(sol_to_plot(13:16,:));
trisurf(hull4,sol_to_plot(13:16,1),sol_to_plot(13:16,2),sol_to_plot(13:16,3),col(4,:)); alpha(.5);
hull5 = convhulln(sol_to_plot(17:20,:));
trisurf(hull5,sol_to_plot(17:20,1),sol_to_plot(17:20,2),sol_to_plot(17:20,3),col(5,:));alpha(.5);
hull6 = convhulln(sol_to_plot(21:24,:));
trisurf(hull6,sol_to_plot(21:24,1),sol_to_plot(21:24,2),sol_to_plot(21:24,3),col(6,:)); alpha(.5);
hold on;


xlabel(['dimension' num2str(dim1)],'Color',[1,0,0],'FontSize',18);  ylabel(['dimension' num2str(dim2)],'Color',[0,1,0],'FontSize',18);
 zlabel(['dimension' num2str(dim3)],'Color',[0,1,0],'FontSize',18);
 
title(['mds on ' subname ' ratings : '   num2str(num_dims_insol) ' dimensional solution'],'FontSize',24);

%find how big the plot needs to be
min_x = min(sol_to_plot(:,dim1));
max_x = max(sol_to_plot(:,dim1));
min_y = min(sol_to_plot(:,dim2));
max_y = max(sol_to_plot(:,dim2));
min_z = min(sol_to_plot(:,dim2));
max_z = min(sol_to_plot(:,dim3));

vals = abs([min_x max_x min_y max_y min_z max_z]);

%largest distance + some wiggle room
val = max(vals) + imsize+1;
%set axes
axis([-val val -val val -val val]);
% colormap([[0:1/255:1]', [0:1/255:1]',[0:1/255:1]'])
set(F,'Position',fig_posn);

hold off;


colormap('default');


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% make a shepard diagram  plotting proximities in the solution against distances in the data 


MakeShepardPlots






