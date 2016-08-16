% converting from rgb to L*a*b*

% 
% want to create an lab space to represent color distributions in.
% 
% lab is a cube which 
% L 0:100  which is black to white
% a and b are opponent axes which go from
% -128:127
% gray values are a=b=0
% 
% the gamut of lab is larger than RGB so we have to find the proportion of Lab that is occupied by RGB and the locations.
% 
% one way to visualize it would simply be to create an rgbcube like we have and then draw the converted points in LAB space


% will need this function
% function [L,a,b] = RGB2Lab(R,G,B)
% RGB2Lab takes red, green, and blue matrices, or a single M x N x 3 image, 
% and returns an image in the CIELAB color space.  RGB values can be
% either between 0 and 1 or between 0 and 255.  Values for L are in the
% range [0,100] while a and b are roughly in the range [-110,110].  The
% output is of type double.
%
% This transform is based on ITU-R Recommendation BT.709 using the D65
% white point reference. The error in transforming RGB -> Lab -> RGB is
% approximately 10^-5.  
%
% See also LAB2RGB.


% so let's make an rgb cube visualization

% need to create equally spaced bins in rgb.  
nbins = 10;
% create the bins.  add 2 to get edges
roundbins = linspace(0,1,nbins+2);
% now find the bin centers.  these will be half the bin size subtracted
% from each bin that is not an edge
roundvec = roundbins(2:end)-roundbins(2)/2;

%   here are all our locations
    bplotmatrix = zeros(length(roundvec)^3,4);
    
%     for each point on red dimension
    for a=1:length(roundvec)
%         for each point on green dimension
        for b=1:length(roundvec)
%             for each point on blue dimension
            for c=1:length(roundvec)
               
%                 fill that point in the 3d matrix with the number of color
%                 matches found with that value

                bplotmatrix(counter,:) = [roundvec(a) roundvec(b) roundvec(c) .1];
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
