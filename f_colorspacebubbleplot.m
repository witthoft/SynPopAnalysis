function f_colorspacebubbleplot(matches, roundvec, dataname,scalefactor)

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
   
    % use downloaded function roundtowardvec
%     roundvec = 0:.25:1; %(1331 locations)
    roundedclrs = roundtowardvec(clrs,roundvec);
    
%     round all colors to the nearest point in the smoothed space
    allroundedclrs = [allroundedclrs; roundedclrs];
    

    counter=1;
    
%     matrix to hold the number of matches at each point in our smoothed
%     space
    bplotmatrix = zeros(length(roundvec)^3,4);
    
%     for each point on red dimension
    for a=1:length(roundvec)
%         for each point on green dimension
        for b=1:length(roundvec)
%             for each point on blue dimension
            for c=1:length(roundvec)
%                 fill that point in the 3d matrix with the number of color
%                 matches found with that value
                bplotmatrix(counter,:) = [roundvec(a) roundvec(b) roundvec(c) ...
                    sum(ismember(roundedclrs, [roundvec(a) roundvec(b) roundvec(c)],...
                    'rows'))];
%                 bplotmatrix is of size length(roundvec)^3 x 4.  the first
%                 3 columns are the 3d coordinates and the fourth column is
%                 the count of matches rounded to the coordinate in that
%                 row.  in this case the bplot matrix is for a single
%                 letter
%             
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
%     this matrix then allows you to sum across all the letters
%   alletterhists will be length(roundvec)^3 * 26 letters  by 4
%   the first 3 columns are the coordinates and the fourth is the count
%     allletterhists = [allletterhists; bplotmatrix];
% let's sum this instead, because I'm not sure how bubbleplot is handling
% repetitions of coordinates

% on the first pass
    if i==1
%         get the coordinates and first set of counts
        allletterhists = bplotmatrix;
    else
%         just add new counts to fourth column
        allletterhists(:,4) = allletterhists(:,4)+bplotmatrix(:,4);
    end
    
    
end




sizes=allletterhists(:,4);

% will need to scale size
% try % of max
% s = bplotmatrix(:,4)/max(bplotmatrix(:,4));

% need to get rid of zeros
p = find(sizes~=0);

figure('name',['distribution across rgb space for ' dataname],'Color',[1 1 1], 'Position',get(0,'ScreenSize'));

subplot(1,2,1);

%     BUBBLEPLOT3(x,y,z,r,c), where c is a rgb-triplet array (in [0,1])
%     with numel(x) rows, plots bubbles with colours specified by c.

% BUBBLEPLOT3(allletterhists(p,1)',allletterhists(p,2)',allletterhists(p,3)
% ',10*(sizes(p)/sum(sizes))',allletterhists(p,1:3));

BUBBLEPLOT3(allletterhists(p,1)',allletterhists(p,2)',allletterhists(p,3)',scalefactor*(sizes(p)/sum(sizes(p)))',allletterhists(p,1:3));

axis equal;
% values range from 1 to 1898
% try using log scaling
% this works but hard to see pattern
% BUBBLEPLOT3(allletterhists(p,1)',allletterhists(p,2)',allletterhists(p,3)',(log(s(p))/(10*max(log(s(p)))))',allletterhists(p,1:3));


% scatter3(roundedclrs(:,1),roundedclrs(:,i,2),roundedclrs(:,i,3))


xlabel('Red');ylabel('Green');zlabel('Blue');
title('data distribution');





% pair it with a uniform distribution
% figure('name',['uniform distribution across rgb space'],'Color',[1 1 1]);

subplot(1,2,2);

% total number of matches
numletterswithmatches = sum(sizes(p));
% number for each point
u = numletterswithmatches/(length(roundvec)^3);

 bplotmatrix(:,4) = u;

 BUBBLEPLOT3(bplotmatrix(:,1)',bplotmatrix(:,2)',bplotmatrix(:,3)',(scalefactor*bplotmatrix(:,4)/numletterswithmatches)',bplotmatrix(:,1:3));
xlabel('Red');ylabel('Green');zlabel('Blue');
set(gca,'Xlim',[0 1],'YLim',[0 1],'ZLim',[0 1]);

axis equal

title('uniform distribution');