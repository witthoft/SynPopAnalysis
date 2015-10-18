
% so the magnet set and the frequency set overlap more than we previously
% thought.  in particular  A,C,H,M,N, & W are all the same
% letter  m fq   high-fq fq no group
% A  red      .81    .69
% C  yellow   .71    .48
% H  orange   .69    .37
% M  red      .58    .37
% N  orange   .75    .31
% W  blue     .62    .29


% could redo the entire analysis on just the subset.

% the selection vector is
nooverlap = [2 4:7 9:12 15:23 25 26];
% no_o_letters = {'B','D','E','F','G','I','J','K','L','O',...
%     'P','Q','R','S','T','U','V','X','Y','Z'}



% it was then updated with some user ids

load EaglemanColoredAlphabets.mat


%  the variable labels gets reused at somepoint... 
letters=labels;
no_o_letters = letters(nooverlap);


% visualize the whole database
fp_visualizeEagDBsubset








% suppose we wanted to something weird, like estimate the expected
% proportions of colors from the labeled database.  that is, what
% percentage of the space is black etc....
load lRGBnathan.mat;





% next we would like to find out how many synesthetes in the database were
% likely to have had the magnet set.  running this script generates a
% number of useful figures

% s_fpEaglemanMatches

s_fpRichAndEaglemanMatches

% another thing is to run this same analysis but comparing subjects to the
% strong group level frequencies in the data.

s_gtEaglemanMatches

s_fqEaglemanMatches
% now let's pull out the subset and graph that.

s_fpRichAndEaglemanMatchessubset


% visualize rich histograms from simulation

% let's do it again as percent
figure('name','percent of times each letter is given a color label rich data', 'Color', [1 1 1]);
for i = 1:26
    
    subplot(6,5,i);
    [counts, bins]=hist(labels.rich(:,i),11)
    counts  = counts/length(labels.rich);
    hBar =bar(bins, counts,'hist');
    ylabel('Number of Subjects');
    %    set(gca,'XTick',[0:length(names)],'XTickLabel',names,'YLim',[0
    %    3000],'FontSize',6);
    set(gca,'YLim',[0 .6],'FontWeight','bold');
    set(hBar,'FaceVertexCData',histcolors);
    title(letters(i));
    box off;
    
end


% this also reminds me that perhaps we should look for subjects who give
% the same color for most letters.  that would be a likely sign of fraud or
% at least unverifiability, and would impact the shuffling in a weird way
% if there were a lot of them.


% we should probably see if we can assign each synesthete to a group and
% generate separate plots for each.  really the point is to analyze the
% subjects left after you take out the frequent template group and the
% magnet group

syntype = zeros(length(dbLabeled),1);


% num of matches needed to be in group
magnetthreshold = 9;
% fqthreshold = 9;

% code for labels variable is in s_fpEaglemanMatches but in this script
% should be in workspace
magmatches = (labelssubset.eagleman == labelssubset.magnets);
% fqmatches = labelssubset.eagleman == labelssubset.fq;
% let's set the fq first
% syntype(find(sum(fqmatches,2)>=fqthreshold))=1;
% then magnet.
syntype(find(sum(magmatches,2)>=magnetthreshold))=2;
% maybe want to find those that are in both groups?

%% plot colors for matches that are close to magnet set


% let's plot each group

% magnets
% matches
figure('name', [ num2str(length(find(syntype==2))) ...
    ' subjects with more than ' num2str(magnetthreshold) ...
    ' matches to letter set'], 'Color', [1 1 1]);

% make a graphical legend 10% as tall as the result
theLegend = rgbsubset.magnets(1:length(find(syntype==2))/10,:,:);
theResult = rgbsubset.eagle(find(syntype==2),:,:);

theStack = [theResult; theLegend];
imagesc(permute(theStack, [1 3 2]))

set(gca, 'XTick', 1:length(nooverlap), 'XTickLabel', no_o_letters, 'YTick', [1 n], 'FontSize', 18)
title([ num2str(length(find(syntype==2))) ...
    ' subjects with more than ' num2str(magnetthreshold) ...
    ' matches to letter set'])


% histogram
% colors come from fp_visualizeEagDB
figure('name','histogram of color names for each letter for magnet syns', 'Color', [1 1 1]);
for i = 1:26
    
    subplot(5,6,i);
    [counts, bins]=hist(dbNumbered(find(syntype==2),i),11)
    counts  = counts/length(find(syntype==2));
    hBar =bar(bins, counts,'hist');
    ylabel('Number of Subjects');
    %    set(gca,'XTick',[0:length(names)],'XTickLabel',names,'YLim',[0
    %    3000],'FontSize',6);
    set(gca,'YLim',[0 1],'FontWeight','bold');
    set(hBar,'FaceVertexCData',histcolors);
    title(letters(i));
    box off;
    
end

% if you wanted to compare it to the earlier figure, you would shuffle the
% data instead of the magnet set and then sort.  the magnet set is just one
% throw instead of n?
% would be better to use the rich, but don't have rgb values, just
% simulated labels

% All the shuffled matches
figure('name', ['top ' num2str(length(find(syntype==2))) ' shuffled matches to letter set'], 'Color', [1 1 1]);

% [Y, ranking] = sort(nummatches.eagleShuffle2magnet, 'descend');
% this is now
[Y, ranking] = sort(nummatches.eagleShuffleByRow, 'descend');

best = ranking(1:length(find(syntype==2)));

% make a graphical legend 10% as tall as the result
theLegend = rgb.magnets(1:round(length(find(syntype==2))/10),:,:);
theResult = rgb.eagleShuffled(best,:,:);

theStack = [theResult; theLegend];
imagesc(permute(theStack, [1 3 2]))

set(gca, 'XTick', 1:26, 'XTickLabel', letters, 'YTick', [1 n], 'FontSize', 18)
title('All shuffled matches to magnet set')







% want to slip a table in here too in the command window
disp('magnet proportions');
disp(['letter  black  white  red    green  yellow blue   brown  purple pink   orange grey']);
for i=1:26
    [counts, bins]=hist(dbNumbered(find(syntype==2),i),11);
    counts  = counts/length(find(syntype==2));
    disp([letters(i) '       ' num2str(counts,'%0.2f   ')]);
    
end


% high fq
figure('name', [ num2str(length(find(syntype==1))) ...
    ' subjects with more than ' num2str(fqthreshold) ...
    ' matches to high fq template'], 'Color', [1 1 1]);

% make a graphical legend 10% as tall as the result
theLegend = rgb.fq(1:length(find(syntype==1))/10,:,:);
theResult = rgb.eagle(find(syntype==1),:,:);

theStack = [theResult; theLegend];
imagesc(permute(theStack, [1 3 2]));

set(gca, 'XTick', 1:26, 'XTickLabel', letters, 'YTick', [1 n], 'FontSize', 18)
title([ num2str(length(find(syntype==1))) ...
    ' subjects with more than ' num2str(magnetthreshold) ...
    ' matches to high fq template']);



% histogram
% colors come from fp_visualizeEagDB
figure('name','histogram of color names for each letter for fq syns', 'Color', [1 1 1]);
for i = 1:26
    
    subplot(5,6,i);
    [counts, bins]=hist(dbNumbered(find(syntype==1),i),11);
    counts  = counts/length(find(syntype==1));
    hBar =bar(bins, counts,'hist');
    ylabel('Number of Subjects');
    %    set(gca,'XTick',[0:length(names)],'XTickLabel',names,'YLim',[0
    %    3000],'FontSize',6);
    set(gca,'YLim',[0 1],'FontWeight','bold');
    set(hBar,'FaceVertexCData',histcolors);
    title(letters(i));
    box off;
    
end

% want to slip a table in here too in the command window
disp('high fq proportions');
disp(['letter  black  white  red    green  yellow blue   brown  purple pink   orange grey']);
for i=1:26
    [counts, bins]=hist(dbNumbered(find(syntype==1),i),11);
    counts  = counts/length(find(syntype==1));
    disp([letters(i) '       ' num2str(counts,'%0.2f   ')]);
    
end


% if you wanted to compare it to the earlier figure, you would shuffle the
% data instead of the magnet set and then sort.  the magnet set is just one
% throw instead of n?  would be good to use rich, but it is only simulated
% labels, not really good to simulate rgb values, though we could

% All the shuffled matches
figure('name', ['top ' num2str(length(find(syntype==1))) ' shuffled matches to high fq template'], 'Color', [1 1 1]);

% [Y, ranking] = sort(nummatches.eagleShuffle2magnet, 'descend');
% is now
[Y, ranking] = sort(nummatches.eagleShuffleByRow, 'descend');
best = ranking(1:length(find(syntype==1)));

% make a graphical legend 10% as tall as the result
theLegend = rgb.fq(1:round(length(find(syntype==1))/10),:,:);
theResult = rgb.eagleShuffled(best,:,:);

theStack = [theResult; theLegend];
imagesc(permute(theStack, [1 3 2]));

set(gca, 'XTick', 1:26, 'XTickLabel', letters, 'YTick', [1 n], 'FontSize', 18);
title(['top ' num2str(length(find(syntype==1))) ' shuffled matches to high fq template']);







% remainder
figure('name', [ num2str(length(find(syntype==0))) ...
    ' subjects with no group'], 'Color', [1 1 1]);

% make a graphical legend 10% as tall as the result
theLegend = rgb.fq(1:length(find(syntype==2))/10,:,:);
theResult = rgb.eagle(find(syntype==0),:,:);

theStack = [theResult];
imagesc(permute(theStack, [1 3 2]));

set(gca, 'XTick', 1:26, 'XTickLabel', letters, 'YTick', [1 n], 'FontSize', 18)
title([ num2str(length(find(syntype==0))) ...
    ' subjects with no group']);



% histogram
% colors come from fp_visualizeEagDB
figure('name','histogram of color names for each letter for no group syns', 'Color', [1 1 1]);
for i = 1:26
    
    subplot(5,6,i);
    [counts, bins]=hist(dbNumbered(find(syntype==0),i),11)
    counts  = counts/length(find(syntype==0));
    hBar =bar(bins, counts,'hist');
    ylabel('Number of Subjects');
    %    set(gca,'XTick',[0:length(names)],'XTickLabel',names,'YLim',[0
    %    3000],'FontSize',6);
    set(gca,'YLim',[0 1],'FontWeight','bold');
    set(hBar,'FaceVertexCData',histcolors);
    title(letters(i));
    box off;
    
end



format short g;

% want to slip a table in here too in the command window
disp('no group proportions');
disp(['letter  black  white  red    green  yellow blue   brown  purple pink   orange grey']);
for i=1:26
    [counts, bins]=hist(dbNumbered(find(syntype==0),i),11);
    counts  = counts/length(find(syntype==0));
    disp([letters(i) '       ' num2str(counts,'%0.2f   ')]);
    
end


% want to slip a table in here too in the command window
disp('magnet proportions');
disp(['letter  black  white  red    green  yellow blue   brown  purple pink   orange grey']);
for i=1:26
    [counts, bins]=hist(dbNumbered(find(syntype==2),i),11);
    counts  = counts/length(find(syntype==2));
    disp([letters(i) '       ' num2str(counts,'%0.2f   ')]);
    
end

% want to slip a table in here too in the command window
disp('high fq proportions');
disp(['letter  black  white  red    green  yellow blue   brown  purple pink   orange grey']);
for i=1:26
    [counts, bins]=hist(dbNumbered(find(syntype==1),i),11);
    counts  = counts/length(find(syntype==1));
    disp([letters(i) '       ' num2str(counts,'%0.2f   ')]);
    
end




% would it be cool to have a 6588x6588 matrix where every entry was the
% number of letters in common between the ith and jth entries?  yes it
% would.


% now for a little fun.  some clustering would be good.  so what we need
% are a representations of the letter set, and then a distance metric, and
% then a clustering method.
% in this case let the representation be the color category
% the distance is the number of letters with different matches (or 26 - num
% matches)

% so dbNumbered is 6588 subjects x 26 letters
% and we need to construct a matrix which is 6588 x 6588 which compares
% each row to every other row

%
% % initialize out matrix
% subjSimMatrix = zeros(6588);
%
% % for each subject
% for i=1:6588
% %     compare every other subject
%     for j=1:6588
% %         get the similarity of those two rows
%         choicesim = length(find(dbNumbered(i,:)-dbNumbered(j,:)==0))
% %         now fill the ith row and the jth column
%         subjSimMatrix(i,j)=choicesim;
%     end
% end
%
%
load subjSimMatrix.mat;

% let's look at the subjects we can't assign to a group
nogroupsim = subjSimMatrix(find(syntype==0),find(syntype==0));
% set diagonal equal to 0
% nogroupsim(logical(eye(size(nogroupsim))))=0;
% % threshold out unrelated syns
% nogroupsim(find(nogroupsim<10))=0;

% so every subject has at least one other 'friend' meaning that they have
% at least 10 matches in common with that person

% let's try some clustering on just this subset
subjDistMatrix = abs(nogroupsim-26);
%
subjSimVect = squareform(subjDistMatrix);
% 
% % let's just look at them all
% for i=6215:length(subjSimMatrix)
% %     make a plot that has all the matches from related subjects and uses
% %     the ith subject as a template
% 
% % All the shuffled matches
% figure('name', ['good matches to subject ' num2str(i)], 'Color', [1 1 1]);
% 
% thrsh = 10;
% % find the index to the subjects in subjSim matrix with more than threshold
% % matches
% frnds = find(subjSimMatrix(i,:)>=thrsh);
% 
% 
% 
% % make a graphical legend 10% as tall as the result
% theLegend = repmat(rgb.eagle(i,:,:),round(length(frnds)/10),1);
% theResult = rgb.eagle(frnds,:,:);
% 
% theStack = [theResult; theLegend];
% imagesc(permute(theStack, [1 3 2]));
% if(length(frnds)>1)
% set(gca, 'XTick', 1:26, 'XTickLabel', letters, 'YTick', [1 length(frnds)], 'FontSize', 18);
% end
% title([num2str(length(frnds)) ' matches for subject ' num2str(i)]);
% 
% saveas(gcf,['tenfriends/subject_' num2str(i) '_10ormore.jpg'],'jpg');
% 
% % input('next');
% close(gcf);
% 
% end

%


%
% % symmetrical matrix in which each entry ranges from 0 to 26, where the
% % entry marks the number of shared matches between subject i and subject j.
%
%
% % let's convert this into the kind of vector output by pdist
% % first turn this into a distance metric
% subjDistMatrix = abs(subjSimMatrix-26);
%
% subjSimVect = squareform(subjDistMatrix);

% 
% if we look at the distribution of birth year for the magnet set we see
% that it fills just the era that the set was produced with a peak in the
% 80s (that is when normalized as percent of the database for that year).
% if we make the same plot for the people who have many of the most
% frequent matches, we would expect that percentage to be uniform across
% time as there is no reason to expect those cultural influences to have
% changed. write a little function that extracts the two groups.
% 
% 
% fix the script that finds the clusters to only pull out the separate
% groups, so set a threshold, and then if there are more than n matches
% make a group that has all of those in it.  then do the remainder.
% 
% 
% need to ask eagleman for the task performance data (test retest
% reliability and speeded decision)
% 
% also classification etc... stuff we had in the first paper
% 
% numbers?
% 
% convert to hsv?
% 
% 
% if you take the people who had the magnet set and look at histograms of
% the letters that don't match, it should look like the rest of the
% population.
% 
% 
% 
% 5
% need true ages of subjects.  right now we just know their age at the time
% of test, but not when they took the test.
% 
% 
% let's try some clustering
% 
%   IDX = KMEANS(X, K) partitions the points in the N-by-P data matrix
%     X into K clusters.  This partition minimizes the sum, over all
%     clusters, of the within-cluster sums of point-to-cluster-centroid
%     distances.  Rows of X correspond to points, columns correspond to
%     variables.  KMEANS returns an N-by-1 vector IDX containing the
%     cluster indices of each point.  By default, KMEANS uses squared
%     Euclidean distances.
% 
% so our matrix is dbNumbered
%   dbNumbered             6588x26
% 
% our problem is choosing a number of clusters which is unknown
% 
% 
% [IDX, C, SUMD, D] = KMEANS(X, K) returns distances from each point
%     to every centroid in the N-by-K matrix D.



% 
% X = dbNumbered(syntype==0);
% K = 2;
% dist = 'cityblock';
% 
% 
% 
% % now we would like to see what kinds of groupings we get in our clusters
% for K=5:25
%     figure('name', ['clusters where K = ' num2str(K)], 'Color', [1 1 1]);
%     [IDX, C, SUMD, D] =kmeans(X,K,'distance',dist);
%     K
%     [nrows ncols] = subplotsize(K);
%     for i=1:K
%         i
%         % for each cluster let's make a figure
%         subplot(nrows,ncols,i);
%         imagesc(p_rgb(find(IDX==i),:,:));
%         ylabel('SUBJECTS');
%         xlabel('LETTERS');
%         % set(gca,'XTick',[1:26],'XTickLabel',labels, 'TickDir','out',
%         % 'YDir','normal');
%         box off;
% 
%     end
% 
% end



% now we can try agglometarive clustering
%  CLUSTERDATA Construct clusters from data.
%     T = CLUSTERDATA(X, CUTOFF) constructs clusters from data X.
%     X is a matrix of size M by N, treated as M observations of N
%     variables.  CUTOFF is a threshold for cutting the hierarchical
%     tree generated by LINKAGE into clusters. When 0 < CUTOFF < 2,
%     clusters are formed when inconsistent values are greater than
%     CUTOFF (see INCONSISTENT). When CUTOFF is an integer and CUTOFF >= 2,
%     then CUTOFF is considered as the maximum number of clusters to
%     keep in the hierarchical tree generated by LINKAGE. The output T is
%     a vector of size M containing a cluster number for each observa
% 
% 
% for clusters = 5:5:50
% 
%     T=clusterdata(X,'distance',dist,'criterion','inconsistent','maxclust',clusters,'linkage','single');
% 
%     figure('name', ['clustering where max is ' num2str(clusters)], 'Color', [1 1 1]);
% 
%     [nrows ncols] = subplotSz(clusters);
%     for i=1:(max(T))
%         % for each cluster let's make a figure
%         subplot(nrows,ncols,i);
%         imagesc(p_rgb(find(T==i),:,:));
%         ylabel('SUBJECTS');
%         xlabel('LETTERS');
%         % set(gca,'XTick',[1:26],'XTickLabel',labels, 'TickDir','out',
%         % 'YDir','normal');
%         box off;
% 
%     end
% 
% end
% 

%
