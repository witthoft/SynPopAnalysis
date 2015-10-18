function    makeIndPlots(matchmatrix,letterstouse)

% function    makeIndPlots(matchmatrix,letterstouse)
% matchmatrix is an m subjects by 26 matrix which has a 1 everywhere the
% subject matches some template and a 0 everywhere else
% letterstouse is an index of numbers between 1 and 26 indiciting what
% subset of letters to do the analysis on
% A B C D E F G H I  J  K L  M   N O   P Q  R  S  T  U  V  W  X  Y  Z
% 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
% should just make that the list of letters.  then you can alter the subset
% code.  duh
% letters is cell array of letters A-Z
letters = ['A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N' 'O' 'P' 'Q' 'R' 'S' 'T' 'U' 'V' 'W' 'X' 'Y' 'Z']';


% % generate P(A&B) matrix
% % so for each column i
for i=1:26
%   check every other column j
    for j=1:26
% we want to know when the ith column and the jth column both matched the
% template.  so essentially where the sum of those columns is 2

        %         for proving indendence we might rather know P(i&j) rather than
        %         P(i|j)  for that we want to divide by the total number of rows
        %         rather than just the match rows
        AandB(i,j)=length(find(matchmatrix(:,i)+...
            matchmatrix(:,j)==2))/length(matchmatrix);
    end
end
% get subset of P(A&B)
AandB = AandB(letterstouse,letterstouse);

% generate base rate matrix
% want the base rate for matches for each letter.
% since you want the rows to predict the column, you can just fill the
% columns with the base rate for that letter to the magnet set
baserate = sum(matchmatrix)/length(matchmatrix);
% now make same size as other matrices
baserate = repmat(baserate,26,1);

% Find P(A)*P(B)
% generate product matrix using base rate of magnet matches and transpose
AxB = baserate(letterstouse,letterstouse).*...
    baserate(letterstouse,letterstouse)';
% unpack uppertriangular parts of P(A&B) and P(A)*_(B) to get these as
% vectors
AxBV=[];
AandBV=[];
for i=1:length(letterstouse)
   AandBV = [AandBV,AandB(i,find((1:length(letterstouse))>i))];
    AxBV = [AxBV,AxB(i,find((1:length(letterstouse))>i))];
end


% get letters for plot
letterpairs = {};

for i=1:length(letterstouse)
    for j=1:length(letterstouse)
        letterpairs{i,j}=[letters(letterstouse(i)) letters(letterstouse(j))];
    end
end

% extract upperT of letter matrix (a bit dumb)
lpairsV=[];
for i=1:length(letterstouse)
    lpairsV = [lpairsV letterpairs(i,find((1:length(letterstouse))>i))];
end


% now we can qualitatively assess the independence
% P(A&B)
figure('Name','Independence Analysis matches to magnet set','Color',[1 1 1]);
subplot(2,2,1);
imagesc(AandB);
h=colorbar;
title('P(A&B)');
set(gca,'CLim',[0 1],'XTick',[1:length(letterstouse)],'YTick',[1:length(letterstouse)],...
    'XTickLabel',letters(letterstouse),'YTickLabel',letters(letterstouse));

% P(A)*P(B)
subplot(2,2,2);
imagesc(AxB);colorbar;title('baserate i x j');
set(gca,'CLim',[0 1],'XTick',[1:length(letterstouse)],'YTick',[1:length(letterstouse)],...
    'XTickLabel',letters(letterstouse),'YTickLabel',letters(letterstouse));
title('P(A)*P(B)');

% histogram of P(A&B)-P(A)*P(B)
subplot(2,2,3);
hist([AandBV-AxBV],[-.1:.025:.1]);
box off;
xlabel('P(A&B)-P(A)*P(B))');
ylabel('number of letter pairs');
title('P(A&B)-P(A)*P(B))');

% scatterplot of P(A&B) vs P(A)P(B)
subplot(2,2,4);
% x =y
plot(0:.01:.5,0:.01:.5,'k--');
axis equal;
hold on;
%
text(AxBV,AandBV,lpairsV);
plot(AxBV,AandBV,'ro');
box off;
set(gca,'XLim',[0 .5],'YLim',[0,.5]);
ylabel('P(A&B)');
xlabel('P(A)*P(B)');
title('P(A&B) vs P(A)*P(B)');
