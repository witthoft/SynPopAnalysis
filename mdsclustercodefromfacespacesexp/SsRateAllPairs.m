% Multidimensional scaling for similarity data collected on 2.08.08.   Very ugly.  Need to load the data (208data.m) which has matrices in the right format.
%For this particular version subjects rated every possible pair of 6 faces
%at 4 views which is 24!/2!(24-2)! = 276 pairs, + every stimulus with
%itself for a total of 300 pairwise distances.   This when finished should
%have complementary plots for both the individual and average data.
%want the cross-correlation matrices
% want the mds solution for each subject and average, as well as subjects
% solutions rotated onto one another
% scree plots to see what the best dimensional explanation might be
% 

clear
% whitebg([1 1 1]);

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


%4.1 tri data
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
% 
% ratings{1}=allsubs; subname = 'average';
% ratings{1}=jonathan; subname = 'jonathan';
% ratings{1}=jin; subname = 'jin';
% ratings{1}=maria; subname='maria';
% ratings{1}=carlo; subname='carlo';
% ratings{1}=charlie; subname='charlie';
% ratings{1}=nina; subname='nina';


% 5.15 tri data

cd /Users/nathan/Experiments/FaceSpaces/ResemblancesWorking/resviewsresults/4.15results
load '415data.mat';
ids = unique(allsubs(:,1));
ratings{1}=allsubs; subname = 'average';
% ratings{1}=erikka; subname = 'erikka';
% ratings{1}=greg; subname = 'greg';
% ratings{1}=steph; subname = 'steph';


% Get id values for the face numbers.  
% ids = unique(allsubst(:,1));
%add second column which is counter (i.e. 1 to number of unique faces)
ids = [ids [1:length(ids)]'];

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
    hist(sims(:),6);
    
    xlabel('rating','FontSize',24);ylabel('number of pairs','FontSize',24);title([subname ' ratings'],'FontSize',24);
    
    
    %start with an empty matrix, and on each pass add another colum
    %corresponding  to a single subject's ratings
%     sims=fliplr(flipud(sims));
    all_sims = [all_sims sims(:)];
    % then Turn similarities into distances.
    dists = 7 - sims;  %rating is from 1 to 7
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
image(triu(sims)); colormap(bone(8)); xlabel('similarities','Color',[1,0,0],'FontSize',18); colorbar;



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
imsize = 0.4;
fig_posn = [0 0 1000 1000];

% Plot overall solution at preset dimensionality (usu 2).
n=n+1;

% F=figure(n);  n=n+1;

%mean_mds_solution is a list of x,y positions
%image is a function which places the center of an matrix  
%image[(x of first face in mds - konst; x of first face in mds +
%konst)], [yo face1 in mds -konst; y of face1 in mds+konst], face1

%so image([left, right],[top, bottom], picture)

% %plot first face
% image([mean_mds_solution(1,1)-imsize; mean_mds_solution(1,1)+imsize],[mean_mds_solution(1,2)-imsize; mean_mds_solution(1,2)+imsize],faces{1});
% hold on;
% 
% %now the rest
% for ctr = 2:length(ids)
%    image([mean_mds_solution(ctr,1)-imsize; mean_mds_solution(ctr,1)+imsize],[mean_mds_solution(ctr,2)-imsize; mean_mds_solution(ctr,2)+imsize],faces{ctr}); 
% end
% 
% text(mean_mds_solution(:,1),mean_mds_solution(:,2),num2str([1:24]'),'Color',[1 1 0], 'FontSize',18);
% 
% xlabel('dimension 1','Color',[1,0,0],'FontSize',18);  ylabel('dimension 2','Color',[0,1,0],'FontSize',18);
% title(['mds on' subname 'data'],'FontSize',24);
% 
% %find how big the plot needs to be
% min_x = min(mean_mds_solution(:,1));
% max_x = max(mean_mds_solution(:,1));
% min_y = min(mean_mds_solution(:,2));
% max_y = max(mean_mds_solution(:,2));
% 
% vals = abs([min_x max_x min_y max_y]);
% 
% %largest distance + some wiggle room
% val = max(vals) + imsize;
% %set axes
% axis([-val val -val val]);
% colormap([[0:1/255:1]', [0:1/255:1]',[0:1/255:1]'])
% set(F,'Position',fig_posn);
% hold off;
% 


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


% plot stress x dimensionality
figure(n);
n=n+1;
plot(strs(:,2),strs(:,1));
xlabel('number of dimensions in solution','FontSize',18);
ylabel('stress','FontSize',18);
title('goodness of fit for solutions of different dimensionality');

% want to plot arbitrary pair of dimensions in a figure
% need to know dimensionality of the solution
% and which two dimensions to plot

num_dims_insol = 2;
dim1 = 1;
dim2 = 2;


%get solution of right dimensionality
sol_to_plot = dim_solutions{num_dims_insol};

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

%now the rest
for ctr = 2:length(ids)
   image([sol_to_plot(ctr,dim1)-imsize; sol_to_plot(ctr,dim1)+imsize],[sol_to_plot(ctr,dim2)-imsize; sol_to_plot(ctr,dim2)+imsize],faces{ctr}); 
end

text(sol_to_plot(1:4,dim1),sol_to_plot(1:4,dim2),num2str([1:4]'),'Color',[1 0 0], 'FontSize',18);
text(sol_to_plot(5:8,dim1),sol_to_plot(5:8,dim2),num2str([5:8]'),'Color',[.5 0 0], 'FontSize',18);
text(sol_to_plot(9:12,dim1),sol_to_plot(9:12,dim2),num2str([9:12]'),'Color',[0 0 1], 'FontSize',18);
text(sol_to_plot(13:16,dim1),sol_to_plot(13:16,dim2),num2str([13:16]'),'Color',[0 0 .5], 'FontSize',18);
text(sol_to_plot(17:20,dim1),sol_to_plot(17:20,dim2),num2str([17:20]'),'Color',[0 1 0], 'FontSize',18);
text(sol_to_plot(21:24,dim1),sol_to_plot(21:24,dim2),num2str([21:24]'),'Color',[0 .5 0], 'FontSize',18);



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

% text(sol_to_plot(:,dim1),sol_to_plot(:,dim2),sol_to_plot(:,dim3),num2str([1:24]'),'Color',[1 0 1], 'FontSize',18);


% text(sol_to_plot(1:4,dim1),sol_to_plot(1:4,dim2),sol_to_plot(1:4,dim3),num2str([1:4]'),'Color',[1 0 0], 'FontSize',18);
% text(sol_to_plot(5:8,dim1),sol_to_plot(5:8,dim2),sol_to_plot(5:8,dim3),num2str([1:4]'),'Color',[.5 0 0], 'FontSize',18);
% text(sol_to_plot(9:12,dim1),sol_to_plot(9:12,dim2),sol_to_plot(9:12,dim3),num2str([1:4]'),'Color',[0 0 1], 'FontSize',18);
% text(sol_to_plot(13:16,dim1),sol_to_plot(13:16,dim2),sol_to_plot(13:16,dim3),num2str([1:4]'),'Color',[0 0 .5], 'FontSize',18);
% text(sol_to_plot(17:20,dim1),sol_to_plot(17:20,dim2),sol_to_plot(17:20,dim3),num2str([1:4]'),'Color',[0 1 0], 'FontSize',18);
% text(sol_to_plot(21:24,dim1),sol_to_plot(21:24,dim2),sol_to_plot(21:24,dim3),num2str([1:4]'),'Color',[0 .5 0], 'FontSize',18);


text(sol_to_plot(1:4,dim1),sol_to_plot(1:4,dim2),sol_to_plot(1:4,dim3),num2str(ids(1:4)'),'Color',[1 0 0], 'FontSize',18);
text(sol_to_plot(5:8,dim1),sol_to_plot(5:8,dim2),sol_to_plot(5:8,dim3),num2str(ids(5:8)'),'Color',[.5 0 0], 'FontSize',18);
text(sol_to_plot(9:12,dim1),sol_to_plot(9:12,dim2),sol_to_plot(9:12,dim3),num2str(ids(9:12)'),'Color',[0 0 1], 'FontSize',18);
text(sol_to_plot(13:16,dim1),sol_to_plot(13:16,dim2),sol_to_plot(13:16,dim3),num2str(ids(13:16)'),'Color',[0 0 .5], 'FontSize',18);
text(sol_to_plot(17:20,dim1),sol_to_plot(17:20,dim2),sol_to_plot(17:20,dim3),num2str(ids(17:20)'),'Color',[0 1 0], 'FontSize',18);
text(sol_to_plot(21:24,dim1),sol_to_plot(21:24,dim2),sol_to_plot(21:24,dim3),num2str(ids(21:24)'),'Color',[0 .5 0], 'FontSize',18);


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
colormap([[0:1/255:1]', [0:1/255:1]',[0:1/255:1]'])
set(F,'Position',fig_posn);
hold off;





cd (thepath);



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

colormap([1 0 0
        .5 0 0
        0 0 1
        0 0 .5
        0 1 0
        0 .5 0]);
caxis([1,6]);

col = [1 1 1 1
        2 2 2 2
        3 3 3 3 
        4 4 4 4
        5 5 5 5
        6 6 6 6];

%now find convex hull for each of the identities
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
%
hold on;

colormap([1 0 0
        .5 0 0
        0 0 1
        0 0 .5
        0 1 0
        0 .5 0]);


sf =1
% 
% text(sol_to_plot(1:4,dim1),sol_to_plot(1:4,dim2),sol_to_plot(1:4,dim3),num2str(ids(1:4)'),'Color',[1 0 0], 'FontSize',18);
% text(sol_to_plot(5:8,dim1),sol_to_plot(5:8,dim2),sol_to_plot(5:8,dim3),num2str(ids(5:8)'),'Color',[.5 0 0], 'FontSize',18);
% text(sol_to_plot(9:12,dim1),sol_to_plot(9:12,dim2),sol_to_plot(9:12,dim3),num2str(ids(9:12)'),'Color',[0 0 1], 'FontSize',18);
% text(sol_to_plot(13:16,dim1),sol_to_plot(13:16,dim2),sol_to_plot(13:16,dim3),num2str(ids(13:16)'),'Color',[0 0 .5], 'FontSize',18);
% text(sol_to_plot(17:20,dim1),sol_to_plot(17:20,dim2),sol_to_plot(17:20,dim3),num2str(ids(17:20)'),'Color',[0 1 0], 'FontSize',18);
% text(sol_to_plot(21:24,dim1),sol_to_plot(21:24,dim2),sol_to_plot(21:24,dim3),num2str(ids(21:24)'),'Color',[0 .5 0], 'FontSize',18);


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


% colormap('default');




% make a shepard diagram  plotting proximities in the solution against distances in the data 















% % Get solutions for each subject.
%  mds_solutions = cell(subjects,1);
% for subject = 1:subjects
%     % Set distance matrix.
%     subject_distances = distances{subject};
%     % Get MDS solution.
%     opts = statset('MaxIter',5000)
% 
%     [mds_dimensions,stress,disparities] = mdscale(subject_distances,soldims,'Criterion','stress','Options',opts);
%     % Store solution.
%     mds_solutions{subject} = mds_dimensions;
% end


% % Go through solutions and find Procrustes transform to fit everything to mean solution as best as possible.
% transformed_solutions = cell(subjects,1);
% for subject = 1:subjects
%     [error,rotated_subject_solution] = procrustes(mean_mds_solution,mds_solutions{subject});
%     transformed_solutions{subject} = rotated_subject_solution;
% end
% 
% 
% 
% % on to plotting stuff
% % For each face, find the bounding co-ordinates.
% bounds = cell(6,1);
% for face = 1:6
%     % Initialize bound to data from the first subject.
%     face_loc = transformed_solutions{1}(face,:);
%     bound = [face_loc(1) face_loc(1) face_loc(2) face_loc(2)];
%     for subject = 1:subjects
%         face_pt = transformed_solutions{subject}(face,:);
%         if face_pt(1) < bound(1)
%             bound(1) = face_pt(1);
%         end
%         if face_pt(1) > bound(2)
%             bound(2) = face_pt(1);
%         end
%         if face_pt(2) < bound(3)
%             bound(3) = face_pt(2);
%         end
%         if face_pt(2) > bound(4)
%             bound(4) = face_pt(2);
%         end
%     end
%     bounds{face} = bound;
% end
% 
% %plot individual solutions in subplots before rotations
% figure(n);n=n+1;
% for i=1:subjects
% %     subplot(1,subjects,i);
%     subplot(2,2,i);
%     axis([-4,4,-4,4]);
% %   plot(mds_solutions{i});
%     text(mds_solutions{i}(:,1),mds_solutions{i}(:,2),num2str([1:24]'),'FontSize',18);
% end
% 


% Plot all transformed solutions.
%this is a plot of all the solutions following procusted transformation
% figure(n); n=n+1;
% set(text, 'FontSize',18)
% text(transformed_solutions{1}(:,1),transformed_solutions{1}(:,2),num2str([1:24]'),'Color','black', 'FontSize',18);
% hold on;
% text(transformed_solutions{2}(:,1),transformed_solutions{2}(:,2),num2str([1:24]'),'Color','blue', 'FontSize',18);
% text(transformed_solutions{3}(:,1),transformed_solutions{3}(:,2),num2str([1:24]'),'Color','red', 'FontSize',18);
% text(transformed_solutions{4}(:,1),transformed_solutions{4}(:,2),num2str([1:24]'),'Color','green', 'FontSize',18);
% % text(transformed_solutions{5}(:,1),transformed_solutions{4}(:,2),num2str([1:6]'),'Color','yellow');
% % text(transformed_solutions{4}(:,1),transformed_solutions{4}(:,2),num2str([1:6]'),'Color',[.50,.50,.50]);
% 
% axis([-4 4 -4 4]);
% % Plot bounding rectangles.
% for face = 1:length(ids)
%     bound = bounds{face};
%     H=fill([bound(1) bound(1) bound(2) bound(2)],[bound(3) bound(4) bound(4) bound(3)],[0 0 0]);
%     set(H,'EdgeColor',[1 1 1]);
%     set(H,'EdgeAlpha',0);
% end
% alpha(0.05);
% hold off;
% 

%same figure but with faces added
%for all views diagram (i.e. for seeing all views for a single subject on
%one plot
% figure(n); n=n+1;
% 
% %load all views
% all_ids = unique(SortAll(:,1));
% %need to sort these from 6 groups of 4 to 4 groups of 6
% 
% all_ids = reshape(all_ids, 4,6)
% all_ids=all_ids'
% all_ids = reshape(all_ids, 24, 1)
% 
% allfaces = cell(1,max(size(all_ids)));
% for ctr = 1:max(size((all_ids)))
%     file = [num2str(all_ids(ctr,1)) '.jpg']
%     allfaces{ctr} = imread(file);
% end
% 
% 
% %combine seperate points into single array
% imsize = 0.3;
% 
% ts = transformed_solutions;
% all_sol=[];
% for i=1:subjects
%     all_sol=[all_sol;ts{i}]
% end
% 
% all_sol = [all_sol(:,1) -1*all_sol(:,2)]
% %plot first face
% image([all_sol(1,1)-imsize; all_sol(1,1)+imsize],[all_sol(1,2)-imsize; all_sol(1,2)+imsize],allfaces{1});
% hold on;
% 
% %now the rest
% for ctr = 2:length(allfaces)
%    image([all_sol(ctr,1)-imsize;all_sol(ctr,1)+imsize],[all_sol(ctr,2)-imsize; all_sol(ctr,2)+imsize],allfaces{ctr}); 
% end
% %find how big the plot needs to be
% min_x = min(all_sol(:,1));
% max_x = max(all_sol(:,1));
% min_y = min(all_sol(:,2));
% max_y = max(all_sol(:,2));
% 
% vals = abs([min_x max_x min_y max_y]);
% axis equal
% 
% hold on;
% set(text, 'FontSize',18)
% text(transformed_solutions{1}(:,1),transformed_solutions{1}(:,2)*-1,num2str([1:6]'),'Color','black', 'FontSize',18);
% text(transformed_solutions{2}(:,1),transformed_solutions{2}(:,2)*-1,num2str([1:6]'),'Color','blue', 'FontSize',18);
% text(transformed_solutions{3}(:,1),transformed_solutions{3}(:,2)*-1,num2str([1:6]'),'Color','red', 'FontSize',18);
% text(transformed_solutions{4}(:,1),transformed_solutions{4}(:,2)*-1,num2str([1:6]'),'Color','green', 'FontSize',18);
% 
% %largest distance + some wiggle room
% val = max(vals);
% %set axes
% axis([-val val -val val]);
% axis equal
% colormap([[0:1/255:1]', [0:1/255:1]',[0:1/255:1]'])
% 
% 
% 



% set(F,'Position',fig_posn);
% hold off;
% axis([-3 3 -3 3]);
% % Plot bounding rectangles.
% for face = 1:length(ids)
%     bound = bounds{face};
%     H=fill([bound(1) bound(1) bound(2) bound(2)],[bound(3) bound(4) bound(4) bound(3)],[0 0 0]);
%     set(H,'EdgeColor',[1 1 1]);
%     set(H,'EdgeAlpha',0);
% end
% alpha(0.05);
% hold off;

%   % Use non-metric scaling to recreate the data in 2D, and make a
%        % Shepard plot of the results.
%        [Y,stress,disparities] = mdscale(dissimilarities,2);
%        distances = pdist(Y);
%        [dum,ord] = sortrows([disparities(:) dissimilarities(:)]);
%        plot(dissimilarities,distances,'bo', ...
%             dissimilarities(ord),disparities(ord),'r.-');
%        xlabel('Dissimilarities'); ylabel('Distances/Disparities')
%        legend({'Distances' 'Disparities'}, 'Location','NorthWest');






% % Plot subject solution.
% figure(subject+5);
% image([mds_dimensions(1,1)-imsize; mds_dimensions(1,1)+imsize],[mds_dimensions(1,2)-imsize; mds_dimensions(1,2)+imsize],faces{1});
% hold on;
% for ctr = 2:20
%     image([mds_dimensions(ctr,1)-imsize; mds_dimensions(ctr,1)+imsize],[mds_dimensions(ctr,2)-imsize; mds_dimensions(ctr,2)+imsize],faces{ctr});
% end
% min_x = min(mds_dimensions(:,1));
% max_x = max(mds_dimensions(:,1));
% min_y = min(mds_dimensions(:,2));
% max_y = max(mds_dimensions(:,2));
% vals = abs([min_x max_x min_y max_y]);
% val = max(vals) + imsize;
% axis([-val val -val val]);
% set(F,'Position',fig_posn);
% title_str = ['Subject ' num2str(subject)];
% title(title_str);
% hold off;

