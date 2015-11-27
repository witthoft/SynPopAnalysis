% get pixelwise similarity between images
% to start with will use sum of squared pixel value differences


% some variables
n=1;


% start by loading images into an m x n x p matrix
% with p images of size m x n

% go to right data set and load images

% 2.08 data
% cd /Users/nathan/Experiments/FaceSpaces/ResemblancesWorking/resviewsresults/2.08.08results/mds
% load 'twoEightSubs.mat'
% ids = unique(Alisha(:,1));


% 2.25 data
% cd /Users/nathan/Experiments/FaceSpaces/ResemblancesWorking/resviewsresults/2.25results/mds
% load '2.25data.mat';
% ids = unique(AllSubs(:,1));


% 3.11.data
% cd /Users/nathan/Experiments/FaceSpaces/ResemblancesWorking/resviewsresults/3.11.results/mds;
% load '311SubData.mat';
% ids = unique(AllSs(:,1));

% 3.24 data
% cd ~/Experiments/FaceSpaces/ResemblancesWorking/resviewsresults/3.24.results/tri/
% load '3.24ratings.mat'
% ids = unique(allsubst(:,1));

%  4.15 data
% cd /Users/nathan/Experiments/FaceSpaces/ResemblancesWorking/resviewsresults/4.15results
% load '415data.mat';
% ids = unique(allsubs(:,1))

% 4.8 rect data
% cd /Users/nathan/Experiments/FaceSpaces/ResemblancesWorking/resviewsresults/4.8results
% load '4.8data.mat';
% ids = unique(allsubs(:,1));

% 8.18 tri data
cd /Users/nathan/Experiments/FaceSpaces/ResemblancesWorking/resviewsresults/8.18results
load '8.18data.mat';
ids = unique(allsubs(:,1))



% ids = testids;


%get shape of subplot
% dimPlot = ceil(sqrt(length(ids)));



%read in image and assign to matrix
%display images in subplot
figure(n);n=n+1;

faces = cell(1,length(ids));

for i = 1:length(ids)
    faces{i}=imread(num2str(ids(i,1)),'jpg');
%     subplot(dimPlot,dimPlot,i);
    subplot(6,4,i);
    imagesc(faces{i}(:,:,1));
    axis equal;
    axis off;
end

% invGray = 1-colormap(gray);
% colormap(invGray);

colormap(gray);

% get size of image
[h,w,d] = size(faces{1});
imagePix = h*w;

% display pixel histograms in subplot
figure(n);n=n+1;
for i = 1:length(ids)
    subplot(6,4,i);
    hist(double(reshape(faces{i}(:,:,1),1,h*w)));
end

% this is a little wacky at the moment since these are jpgs and therefor
% are m x n x 3 images.   however, since they are greyscale, all 3 layers
% have the same value. another thing is that the jpgs get read in as uint8
% which means they only have a range from 0-255 (so 239-255=0)

% lets change the representation of the faces to make them n cases x p
% variables and then convert them to a matrix
% matrix is big due to number of pixels, may want to do some averaging to
% something like 400 or so.  also convert those unit8s to doubles.
facematrix=[];
for i=1:length(ids)
    facematrix=[facematrix;double(reshape(faces{i}(:,:,1),1,h*w))];
end

% unjustified hack:  going to switch range from 0-255 to 0-7
facematrix = facematrix*7/255;


%compute sum of squared pixel differences
for i = 1:length(ids)
    for j = 1:length(ids)

        PixDist(i,j)=sum((facematrix(i,:)-facematrix(j,:)).^2)/imagePix;
        
%         PixDist(i,j) = sqrt(sum(sum((normfaces{i}(:,:,1)-normfaces{j}(:,:,1)).^2)))/imagePix;
    end
end



figure(n);n=n+1;
% imagesc(triu(PixDist));colormap(hot);
imagesc(PixDist);colormap(hot);
 colorbar;  
hold off;


figure(n);n=n+1;
 plot(PixDist(find(PixDist(25:48))))
 
% now need to write an output that is useful for mds and hc tools I wrote
% these should be of the form face1 face 2 similarity.  should fix with
% flags in the future, but for now just want to see what happens.

% convert distances to similarities
% as a hack should scale these to range from 1-7 just like face similarity
% data.  this means I can get the faces to shown on the plots

PixSims = abs(PixDist-max(max(PixDist)))

PixSimData=[];
for i=1:length(ids) %sample row
    for j=i:length(ids) %sample uppertriangular columns
        PixSimData=[PixSimData; ids(i) ids(j) PixSims(i,j)];
    end
end
