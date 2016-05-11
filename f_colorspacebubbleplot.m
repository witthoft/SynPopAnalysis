function f_colorspacebubbleplot(matches, roundvec, dataname,scalefactor,circparam)

% makes a 3d bubbplot of showing the proportion of matches for synesthesia
% data at each point in rgb space.  does it for all matches and each letter
% matches are nx26x3 matrix of matches
% binning bins the rgb space into managable size
% dataname puts a name on the figure and could be used for saving




% this is always handy for labelling not sure if used
letters = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N' 'O' 'P' 'Q' 'R' 'S' 'T' 'U' 'V' 'W' 'X' 'Y' 'Z'};






% counts of letters in each rounded 3d bin across all letters
% this is what will be passed to bubbleplot once we fill it
allletterhists = [];

% rounded rgb values for color matches.  basically rounded into the nearest
% bin
allroundedclrs = [];

% counter to keep track of where we are in the loop
clrctr = 1;

% for each letter in matches get the number of matches for 
for i=1:size(matches,2)
%     matches for all subjects for letter i
    clrs = squeeze(matches(:,i,:));
   
    % use downloaded function roundtowardvec
%     rounds each rgb coordinate to nearest location in roundvec
% 
    roundedclrs = roundtowardvec(clrs,roundvec);

    
%     round all colors to the nearest point in the smoothed space
% basically concatenate rounded colors across letters ino a single rgb
% matrix  

    allroundedclrs = [allroundedclrs; roundedclrs];
    
    
    
    
    
%     
    
    counter=1;
    
%     matrix to hold the number of matches at each point in our smoothed
%     space  seems like this could have been done just once but ....
%   bplotmatrix is xyz (rgb) and then the count of matches for each
%   coordinate
%   set that matrix to zeros
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
   
    

%     for reasons which are puzzling to me we are building up
%     allletterhists by summing the bubbleplotmatrix across each letter.  I
%     think we could have once?  anyway

% on the first pass
    if i==1
%         get the coordinates and first set of counts
        allletterhists = bplotmatrix;
    else
%         just add new counts to fourth column
        allletterhists(:,4) = allletterhists(:,4)+bplotmatrix(:,4);
    end
    
    
end

% at the end of this we have allletterhists which is r g b count x roundvec
% length


% sizes is just the count for each rgb coordinate.  not sure why I pulled
% it out and made it its own variable.  maybe shorter to write?
sizes=allletterhists(:,4);

% need to get rid of zeros
p = find(sizes~=0);

% there are two scaling problems

% the first is setting the size of the bubbles so that for whatever binning
% of the space we use (that is how finely we chop up RGB) our bubbles will
% not be so big they overlap with all the others or so small we can't see
% them.  this is handeled by the variable scale factor which is passed to
% this function

% the second is scaling the bubbles so that the size differences relative
% to say uniform can be reasonably understood.  for example if the matches
% were random uniform distributed, than each bubble would have let's say a
% circumference of 1 (bubble plot sizes the bubbles by treating the
% count as setting the radius).  but I think what the eye sees is
% the area of the circle.  suppose we want to make it so the area goes up
% linearly with the count.  
% the easiest way is to say our counts are areas and then find the
% corresponding circumference for each area and pass that to bubbleplot

% area = C^2/(4pi) so 
% sqrt(4piA)=C

% probably should turn this into a switch

switch(circparam)
    case('radius')
        %         do nothing
    case('area')
        % linear with area
        sizes = sqrt(sizes*pi);
    case('diameter')
        % linear with diameter
        sizes = sizes*2
end
% allletterhists(:,4) = sqrt(allletterhists(:,4)*4*pi);
% fucker moved this to sizes

% default is linear with circumference



% linear with volume







figure('name',['distribution across rgb space for ' dataname],'Color',[1 1 1], 'Position',get(0,'ScreenSize'));

subplot(1,2,1);

%     BUBBLEPLOT3(x,y,z,r,c), where c is a rgb-triplet array (in [0,1])
%     with numel(x) rows, plots bubbles with colours specified by c.

% BUBBLEPLOT3(allletterhists(p,1)',allletterhists(p,2)',allletterhists(p,3)
% ',10*(sizes(p)/sum(sizes))',allletterhists(p,1:3));

bubbleplot3(allletterhists(p,1)',allletterhists(p,2)',allletterhists(p,3)',scalefactor*(sizes(p)/sum(sizes(p)))',allletterhists(p,1:3));

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

bubbleplot3(bplotmatrix(:,1)',bplotmatrix(:,2)',bplotmatrix(:,3)',(scalefactor*bplotmatrix(:,4)/numletterswithmatches)',bplotmatrix(:,1:3));
xlabel('Red');ylabel('Green');zlabel('Blue');
set(gca,'Xlim',[0 1],'YLim',[0 1],'ZLim',[0 1]);

axis equal

title('uniform distribution');