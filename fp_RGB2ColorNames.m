% script for dividing up rgb space into color categories.


% so the rgb space is all possible combinations of 3 256 element vectors.
% so 256^3 =  16,777,216 possible colors which is a lot


% we would like to assign every element in that space to an arbitrary
% number of color categories.  one way would be to write down the name for
% every color which is insane.  it would actually be faster to hand label
% every match in the eagleman database!



% berlin and kay color categories: black, blue, brown, grey, green, orange,
% pink, purple, red, white, and yellow


% so want to assign every point in the rgb cube one of these 11 color names

% what is an easy way to do that?

% one thing is that we only need to find the boundaries.  so for example,
% if we just looked at the 255 slices where the bottom slice is b = 0 (that
% is where we make r->x, g->y, z->b) then you could imagine a graphical
% interface where you label the boundaries with a color name
% (or better a pair of color names representing the boundary?) and then the
% interior is automatically completed.

% make rgb matrix

rgb = [0:255;0:255;0:255]';

% make a matrix for color categories (this is a cube);

colornames = zeros(256,256,256);

% each rgb value is then an index to colornames.


% and then we need a readout table


% look at an rgb slice at %bv blue value
bv=.1;
figure;
r = repmat([0:255]',1,256);
r=r/255;
g= r';
b=bv*ones(256,256);

slice = cat(3, r, g, b);
imagesc(slice);



figure;
% try down sampling sample look at an rgb slice
for i=0:10
bv=i*.1;
subplot(2,6,i+1);
r = repmat([0:32:255]',1,8);
r=r/256;
g= r';
b=bv*ones(8,8);

slice = cat(3, r, g, b);
imagesc(slice);

end

% here we have 11 blue values and 16 x 16 resolution in the red/green
% dimension which is 2816 color values


% could also do 8*8*10 which would be 640 trials

% I should probably just do both to see if one turns out better than the
% other



% so the easiest thing is just to do some low resolution naming.



% so we want a program where we can submit a mxnxp (r g b) matrix
% we also want to submit a list of color choices to assign to each entry
% and their corresponding numbers

%  it should show all the entries in random order and allow the user to
%  supply a label for each entry

% output is the color list and labeled matrix and a version of that matrix
% which is scaled to the full color space, though that could be a separate
% function.  it may not be necessary to scale up, as you could down sample
% the color vectors submitted for classification



% 





