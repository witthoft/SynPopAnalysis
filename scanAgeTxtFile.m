% convert text file with ages and dobs to matlab matrix

% the original data file is a bit hard to parse so I pulled out the usable
% subjects using excel.  there were various problems, such as subjects giving
% no birthdate and no age, giving an age but all having done the test in
% 1969, giving impossible ages or birth dates, writing numbers as text (for
% example fourteen for 14, or giving inexact ages such as about 50.  
% open the text file
fID = fopen('AgeUsableSubs.txt');

% get it into matlab
AgeDataTable=textscan(fID,'%d %d %d %d');

% close the file
fclose(fID);

% the columns are userid  batteryid  yearborn and age in 2014 might be off
% by 1 year.  many of these subjects gave a birth year, the rest I inferred
% using their age at test and the date of the text


%convert to numbers
AgeMatrix = [AgeDataTable{1} AgeDataTable{2} AgeDataTable{3} AgeDataTable{4}];
save('AgeMatrix.mat','AgeMatrix');

% find the subjects that we are actually using



% get those into the same order as the rgb colors



% double check using the data we already have