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
% %
% % % let's look at the subjects we can't assign to a group
% % nogroupsim = subjSimMatrix(find(syntype==0),find(syntype==0));
% % set diagonal equal to 0
% % nogroupsim(logical(eye(size(nogroupsim))))=0;
% % % threshold out unrelated syns
% % nogroupsim(find(nogroupsim<10))=0;
% %
% % % so every subject has at least one other 'friend' meaning that they have
% % % at least 10 matches in common with that person
% %
% % % let's try some clustering on just this subset
% % subjDistMatrix = abs(nogroupsim-26);
% % %
% % subjSimVect = squareform(subjDistMatrix);
% %
% % % let's just look at them all
% % for i=6215:length(subjSimMatrix)
% % %     make a plot that has all the matches from related subjects and uses
% % %     the ith subject as a template
% %
% % % All the shuffled matches
% % figure('name', ['good matches to subject ' num2str(i)], 'Color', [1 1 1]);
% %
% % thrsh = 10;
% % % find the index to the subjects in subjSim matrix with more than threshold
% % % matches
% % frnds = find(subjSimMatrix(i,:)>=thrsh);
% %
% %
% %
% % % make a graphical legend 10% as tall as the result
% % theLegend = repmat(rgb.eagle(i,:,:),round(length(frnds)/10),1);
% % theResult = rgb.eagle(frnds,:,:);
% %
% % theStack = [theResult; theLegend];
% % imagesc(permute(theStack, [1 3 2]));
% % if(length(frnds)>1)
% % set(gca, 'XTick', 1:26, 'XTickLabel', letters, 'YTick', [1 length(frnds)], 'FontSize', 18);
% % end
% % title([num2str(length(frnds)) ' matches for subject ' num2str(i)]);
% %
% % saveas(gcf,['tenfriends/subject_' num2str(i) '_10ormore.jpg'],'jpg');
% %
% % % input('next');
% % close(gcf);
% %
% % end
% %
% % %
% %
%
% %
% % % symmetrical matrix in which each entry ranges from 0 to 26, where the
% % % entry marks the number of shared matches between subject i and subject j.
%
%
% % let's convert this into the kind of vector output by pdist
% % first turn this into a distance metric
subjDistMatrix = abs(subjSimMatrix-26);
% %
subjSimVect = squareform(subjDistMatrix);
%
% %
% % if we look at the distribution of birth year for the magnet set we see
% % that it fills just the era that the set was produced with a peak in the
% % 80s (that is when normalized as percent of the database for that year).
% % if we make the same plot for the people who have many of the most
% % frequent matches, we would expect that percentage to be uniform across
% % time as there is no reason to expect those cultural influences to have
% % changed. write a little function that extracts the two groups.
% %
% %
% % fix the script that finds the clusters to only pull out the separate
% % groups, so set a threshold, and then if there are more than n matches
% % make a group that has all of those in it.  then do the remainder.
% %
% %
% % need to ask eagleman for the task performance data (test retest
% % reliability and speeded decision)
% %
% % also classification etc... stuff we had in the first paper
% %
% % numbers?
% %
% % convert to hsv?
% %
% %
% % if you take the people who had the magnet set and look at histograms of
% % the letters that don't match, it should look like the rest of the
% % population.
% %
% %
% %
% % 5
% % need true ages of subjects.  right now we just know their age at the time
% % of test, but not when they took the test.
%
% %
% % let's try some clustering
% %
% %   IDX = KMEANS(X, K) partitions the points in the N-by-P data matrix
% %     X into K clusters.  This partition minimizes the sum, over all
% %     clusters, of the within-cluster sums of point-to-cluster-centroid
% %     distances.  Rows of X correspond to points, columns correspond to
% %     variables.  KMEANS returns an N-by-1 vector IDX containing the
% %     cluster indices of each point.  By default, KMEANS uses squared
% %     Euclidean distances.
% %
% % so our matrix is dbNumbered
% %   dbNumbered             6588x26
% %
% % our problem is choosing a number of clusters which is unknown
% %
% %
% % [IDX, C, SUMD, D] = KMEANS(X, K) returns distances from each point
% %     to every centroid in the N-by-K matrix D.
%
%
%
%
X = dbNumbered;
K = 2;
dist = 'cityblock';
%
%
for i=1:10
    % % now we would like to see what kinds of groupings we get in our clusters
    % % for K=10:2:40
    for K=9
        %     make a figure
        figure('name', ['clusters where K = ' num2str(K)], 'Color', [1 1 1]);
        %     compute k means solution of size k
        [IDX, C, SUMD, D] =kmeans(X,K,'distance',dist);
        %     IDX is number of cluster each subject belongs to and will range from
        %     1:k
        %    C is the center of each cluster in 26 d space (so C is Kx26)
        %     SUMD is the sum of all the distances in a cluster to its center
        %    D is the distance of every point to every center and is 6588 x K
        
        %     compute a fit metric for this k
        %    cna use the average mahalanobis distance with the following algorithm
        %    X is your data
        %    C are your cluster Centers
        %    going to use Y=p/2 where p is dimensionality of data
        %    the justification for this would take me a while
        
        %    the algorighm is
        %    1.  given X generate K clusters
        %    2.  compute a distortion metric d which is
        %     this turns out to be hard!  need to spend more time.
        %  for now let's try using D
        
        
        %   set up subplot of kmeans clusters
        [nrows ncols] = subplotsize(K);
        %     for each cluster
        for i=1:K
            % make a figure showing the letter color matches of its members
            subplot(nrows,ncols,i);
            imagesc(p_rgb(find(IDX==i),:,:));
            ylabel('SUBJECTS');
            xlabel('LETTERS');
            set(gca,'XTick',[1:26],'XTickLabel',letters );
            % 'YDir','normal');
            box off;
            
        end
        
    end
    
end
Kdists = zeros(1,20);

% let's try to find an optimal K
for K=1:length(Kdists)
    
    %     compute k means solution of size k
    [IDX, C, SUMD, D] =kmeans(X,K,'distance',dist);
    %     IDX is number of cluster each subject belongs to and will range from
    %     1:k
    
    Kdists(K)=mean(SUMD);
    %    C is the center of each cluster in 26 d space (so C is Kx26)
    %     SUMD is the sum of all the distances in a cluster to its center
    %    D is the distance of every point to every center and is 6588 x K
    
    %     compute a fit metric for this k
    %    cna use the average mahalanobis distance with the following algorithm
    %    X is your data
    %    C are your cluster Centers
    %    going to use Y=p/2 where p is dimensionality of data
    %    the justification for this would take me a while
    
    %    the algorighm is
    %    1.  given X generate K clusters
    %    2.  compute a distortion metric d which is
    %     this turns out to be hard!  need to spend more time.
    %  for now let's try using D
    
    
    
    
end

% make a plot of average distance from cluster center

figure('Name','average distance from center as a function of K','Color',[1 1 1]);
plot(4:length(Kdists),abs(diff(Kdists(3:end))),'ro');
box off;
xlabel('number of clusters'); ylabel('decrease in average distance with additional cluster');




% one problem with k-means is that it expects clusters to be of similar
% sizes.  given that we don't expect that, it would be better to use a
% mixture of gaussians model to find our clusters.


% what the hell is a mixture of gaussians model? try to figure that out...


% there is a thing called the density estimation problem.  you have many
% points in a D-dimensional space.  the points may be grouped in clusters.

% so you have n points in D dimensions.
% those points were generated by some unknown number of probability
% distributions RD from family F.

% the problem is to find the f(x) element of F that generated your data


% that seems hard!


% Assume that the functions in F are mixtures of gaussians. so each
% function is N-dimensional with parameters specifying the distribution on
% each dimension.  Gaussians sum to 1 so these are all probability
% distributions.

% what you are doing I think is fitting all of these gaussians at once to
% the data.  each one has its own mean and covariance.  and each one has a
% weight, so that all the gaussians together sum to 1.
% I'm guessing that each gaussian captures a cluster, and that clustering
% is defined probabilistically, because each actual data point has a
% probaability of being in each cluster (soft clustering).

% in matlab we fit by doing

% obj = gmdistribution.fit(X,k,...,param1,val1,param2,val2,...)

% where X is our data and k is the number of gaussians.  like k-means we
% fit separately for each proposed number of sources.


%
%
% % some tricks about this
% % fit model
% o = gmdistribution.fit(X,K,'Regularize',.1)
%
% % can assign each subject to one of k clusters using
%
%
% for K=2:15
% %     make a figure
%     figure('name', ['mixtures of gaussians: clusters where K = ' num2str(K)], 'Color', [1 1 1]);
% %     compute k means solution of size k
% %
% %     o = gmdistribution.fit(X,K,'Regularize',.00001)
%     o = gmdistribution.fit(X,K,'SharedCov',true)
%
%     IDX = cluster(o,X);
%
% %   set up subplot of kmeans clusters
%     [nrows ncols] = subplotsize(K);
% %     for each cluster
%     for i=1:K
%         % make a figure showing the letter color matches of its members
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
%
%
%




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