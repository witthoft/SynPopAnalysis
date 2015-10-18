% want to make a csv file for dan like this
%
% any chance you could give it to me as a CSV where each line / row looks like this?
%
% subject_id,letter,rval,gval,bval


cd '/Volumes/biac4-kgs/Projects/synaesthesia/Eagleman'


load RGBEagle.mat

% whos
%   rgbEagle      6588x3x26            4110912  double

% new matrix to make
EagleRGB =[];

rowcount = 1;
% for each subject
for s=1:length(rgbEagle)
    %  for each letter
    for l=1:26
%     write that row of the matrix
    EagleRGB(rowcount,:) = [s l rgbEagle(s,:,l)];
%     move to next row
    rowcount=rowcount+1;
    end

end

csvwrite('RGBEaglmanDB.csv',EagleRGB);


% check on other data

for l=1:26
[n c]=hist(labels.eagleman(:,l),0:11);
letters(l)
flipud(sort(n)')
end
