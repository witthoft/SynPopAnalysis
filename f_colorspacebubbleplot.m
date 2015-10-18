function f_colorspacebubbleplot(matches, roundvec, dataname)

% makes a 3d bubbplot of showing the proportion of matches for synesthesia
% data at each point in rgb space.  does it for all matches and each letter
% matches are nx26x3 matrix of matches
% binning bins the rgb space into managable size
% dataname puts a name on the figure and could be used for saving


letters = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N' 'O' 'P' 'Q' 'R' 'S' 'T' 'U' 'V' 'W' 'X' 'Y' 'Z'};





% counts of letters in each rounded 3d bin across all letters
allletterhists = [];

% rounded rgb values across all letters
allroundedclrs = [];
clrctr = 1;

% for each letter
for i=1:26
    clrs = squeeze(matches(:,i,:));
    %  all rgb values are between 0 and 1
    % we could round to the nearest .1
    % use downloaded function roundtowardvec
%     roundvec = 0:.25:1; %(1331 locations)
    roundedclrs = roundtowardvec(clrs,roundvec);
    
    allroundedclrs = [allroundedclrs; roundedclrs];
    
    % want to be able to scale our plots (either alpha or markersize) using the
    % number of datapoints in each bin.  3d hist code in matlab looks like a
    % pain so we could do it in a loop?
    
    counter=1;
    bplotmatrix = zeros(length(roundvec)^3,4);
    
    for a=1:length(roundvec)
        for b=1:length(roundvec)
            for c=1:length(roundvec)
                bplotmatrix(counter,:) = [roundvec(a) roundvec(b) roundvec(c) ...
                    sum(ismember(roundedclrs, [roundvec(a) roundvec(b) roundvec(c)],...
                    'rows'))];
                counter=counter+1;
            end
        end
    end
   
    
%     can add letterwise plots back in later
%     figure('name',['distribution across rgb space for ' dataname ' ' letters{i}],'Color',[1 1 1]);
%    
%     sizes=bplotmatrix(:,4);
%     
%  
%     % need to get rid of zeros
%     p = find(sizes~=0);
%     
%     %     BUBBLEPLOT3(x,y,z,r,c), where c is a rgb-triplet array (in [0,1])
%     %     with numel(x) rows, plots bubbles with colours specified by c.
%     
%     BUBBLEPLOT3(bplotmatrix(p,1)',bplotmatrix(p,2)',bplotmatrix(p,3)',(sizes(p)/sum(sizes))',bplotmatrix(p,1:3));
%     % scatter3(roundedclrs(:,1),roundedclrs(:,i,2),roundedclrs(:,i,3))
%     
%     
%     xlabel('Red');ylabel('Green');zlabel('Blue');
%     
    %     collect everything into a giant matrix
    
    allletterhists = [allletterhists; bplotmatrix];
    
    
end




sizes=allletterhists(:,4);

% will need to scale size
% try % of max
% s = bplotmatrix(:,4)/max(bplotmatrix(:,4));

% need to get rid of zeros
p = find(sizes~=0);

figure('name',['distribution across rgb space for ' dataname],'Color',[1 1 1]);


%     BUBBLEPLOT3(x,y,z,r,c), where c is a rgb-triplet array (in [0,1])
%     with numel(x) rows, plots bubbles with colours specified by c.

BUBBLEPLOT3(allletterhists(p,1)',allletterhists(p,2)',allletterhists(p,3)',10*(sizes(p)/sum(sizes))',allletterhists(p,1:3));

% values range from 1 to 1898
% try using log scaling
% this works but hard to see pattern
% BUBBLEPLOT3(allletterhists(p,1)',allletterhists(p,2)',allletterhists(p,3)',(log(s(p))/(10*max(log(s(p)))))',allletterhists(p,1:3));


% scatter3(roundedclrs(:,1),roundedclrs(:,i,2),roundedclrs(:,i,3))


xlabel('Red');ylabel('Green');zlabel('Blue');






% pair it with a uniform distribution
figure('name',['uniform distribution across rgb space'],'Color',[1 1 1]);

% total number of matches
numletterswithmatches = sum(sizes(p));
% number for each point
u = numletterswithmatches/(length(roundvec)^3);

 bplotmatrix(:,4) = u;

 BUBBLEPLOT3(bplotmatrix(:,1)',bplotmatrix(:,2)',bplotmatrix(:,3)',(10*bplotmatrix(:,4)/numletterswithmatches)',bplotmatrix(:,1:3));
xlabel('Red');ylabel('Green');zlabel('Blue');


