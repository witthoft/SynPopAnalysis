% function [labeledRGB colornames] = LabelRGBIndices(matrixSize,
% colornames, subjectname)

% [lRGB colornames] = LabelRGBIndices([8,8,11])
% did this.  now can load the labeled matrix

load lRGBnathan.mat

% labeld matrix is called labeldRGB and is 9x9x12

% now need a script that plots the colors and my choices so I can check how
% it went

% here are our original color indices
matrixSize = [8 8 11];

% red
r = [0:256/matrixSize(1):256];
% green
g = [0:256/matrixSize(2):256];
% blue
b = [0:256/matrixSize(3):256];


% so for each blue level want to show on one side the original rgb color
% and then on the other the matched choice.

% here are some good examples for each color

% need a series of rects for our match buttons
black = [0 0 0];
white = [255 255 255];
red = [255 0 0];
green = [0 255 0];
yellow = [255 255 0];
blue = [0 0 255];
brown = [139 69 19];
purple = [148 0 211];
pink = [255 20 147];
orange = [255 165 0];
gray = [150 150 150];

% should stack these into a single matrix but for reading clarity will
% leave separate even though code is longer
%     so we will have 11 color choices
%     0 = black
%     1 = white
%     2 = red
%     3 = green
%     4 = yellow
%     5 = blue
%     6 = brown
%     7 = purple
%     8 = pink
%     9 = orange
%    10 = grey


%  so step through our labeled matrix and then make an rgb image for each
%  slice

% actually let's combine to make the indexing easier
colors = [black;white;red;green;yellow;blue;brown;purple;pink;orange;gray]/255;

lRGB=labeledRGB+1;
figure;
% make a series of blue slices
for i=1:size(labeledRGB,3)
    
    set(gcf,'colormap',colors);
    subplot(4,3,i);
    axis equal
    image(lRGB(:,:,i));
    axis equal
    axis tight
    
end

% ok now lets generate the actual rgb colors

% red
r = [0:256/8:256]';
% green
g = [0:256/8:256];
% blue
b = [0:256/11:256];

% make the red and green matrices which are orthogonal
re=repmat(r,1,length(g))/256;
gr=re'
bl=ones(size(re))/256;

figure;
for i=1:length(b)
    subplot(4,3,i);
    slice = cat(3, re, gr, (b(i))*bl);
    imagesc(slice);
    axis equal;
    axis tight;
end




% can compare and see that things look pretty good.
% so what do we want to do with this?
% we would like to take a list of rgb values from say 1-26 for the alphabet
% interpolate to our rgb locations and then return 26 color names.   seems
% simple

% but of course it is more confusing than I thought.  one way to think
% about it is to ask what is the nearest location in the downsampled cube
% and then just read out the color value there.  but another way is to
% consider the space to be r, g, b, and then find the nearest point in that
% space. its the same but the indexing is slightly different.




% looks like we can use dsearchn.  which I should figure out how that works

% first we have to produce all possible rgb values in our matrix

rgbgrid = [];
% dumb but effective?
for i=1:length(r)
    for j=1:length(g)
        for k=1:length(b)
            rgbgrid = [rgbgrid; r(i) g(j) b(k)];
        end
    end
end

figure;

for i=1:100
    % easy.  so let's say we have a point
    p=255*rand(1,3);
    %  let's look at what color it is:
    
    %     figure('Color',p/256);
    subplot(10,10,i);
    %     this is retarded; should just be generated in the parentheses.  why
    %     even comment instead of doing it?
    c = p/256;
    imagesc(reshape(c,1,1,3));
    
    axis off;
    % this gives us the row in our colorlist
    newp = dsearchn(rgbgrid,p);
    
    % this gives us the entries in the row
    intrgb = rgbgrid(newp,:);
    
    % this gives us the location in the labeled rgb cube
    
    labelloc = [find(r==intrgb(1)) find(g==intrgb(2)) find(b==intrgb(3))];
    
    % and now the color that goes with that location
    
    colornames = {'black', 'white', 'red', 'green','yellow','blue', 'brown',...
        'purple', 'pink', 'orange', 'gray'};
    
    % since I was not planning ahead, my label numberings go from 0-10 instead
    % of 1-11 which means I have to add 1 to pick the right number out of the
    % list
    title(colornames(labeledRGB(labelloc(1),labelloc(2),labelloc(3))+1));
end

% running this shows maybe 5-15% error rate in the labelling.  most of the
% mixes are green/blue  green/yellow  and pink/red.  anyway, need to do
% better but this works in principle


% let's test our labelled matrix on our synesthesia data.
load synstruct.mat
% loads a struct called syn which has all the data from our psych science
% paper
% this was the coding I used for the first paper. which of course was
% totally different.
%  1=red, 2=orange, 3=yellow, 4=green, 5=blue, 6=purple,
%  7=brown, 8=pink  9=black, 10=white, 11=gray,

pscicolornames = {'red' 'orange' 'yellow' 'green' 'blue' 'purple' 'brown' ...
    'pink' 'black' 'white' 'gray'};






for i=1:length(syn)
    %     for each subject make a figure
    figure('Name',syn(i).name);
    for j =1:length(syn(i).rgb)
        %         then plot each of their matches with the label I gave on the top
        %         and the label given from the matrix on the bottom
%         scale our color for ease
        c=syn(i).rgb(j,:)/255;
        subplot(5,6,j);
        imagesc(reshape(c,1,1,3));
        %        label I gave
        title(pscicolornames(syn(i).cats(j)));
        %         label from matrix
        % this gives us the row in our colorlist
        newp = dsearchn(rgbgrid,c*255);
        
        % this gives us the entries in the row
        intrgb = rgbgrid(newp,:);
        
        % this gives us the location in the labeled rgb cube
         labelloc = [find(r==intrgb(1)) find(g==intrgb(2)) find(b==intrgb(3))];
        
        xlabel(colornames(labeledRGB(labelloc(1),labelloc(2),labelloc(3))+1));
        
    end
    
    
end



% this works pretty good.  in fact it looks like there are a couple of data
% entry errors in the syn struct for the category labels as the automatic
% labelling is correct but the hand label is incorrect.  not enough to
% change anything.  there are some disagreements about brown and orange,
% where the color is probably best labelled as brown but is so close to
% orange that it fits into the letterset pattern nicely.  that is simply a
% drawback of this method I suppose.  some similarity relations will be
% lost.

% I mean you could start representing the color matches in terms of
% direction and magnitude in the color space.  so each letter would be
% represented by some polar coordinate relative to the preceding one.
% probably not worth doing, but the question of how to represent the data
% in the best way to find similar sets is an important one.  

