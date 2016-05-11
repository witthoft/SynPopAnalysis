% want to use bigram frequencies to predict color distances for
% synesthetes.  will us the bigrams that come from google
% I got the raw count of word occurences from
% http://norvig.com/mayzner.html
% which he got from here:
% http://storage.googleapis.com/books/ngrams/books/datasetsv2.html


% so this script should
% 1. load up the n most frequent words
% 2. for every letter pair find how many times it occurs
% 3. convert those counts to frequencies


% number of words we want to consider
% this may not matter, but since the vocabulary of a 6 year old is smaller
% than an adult etc...
% my guess is that above some number the frequencies don't change that much

% let's start with 10,0000
nwords = 5000;


% load the data

fname = 'WORDCOUNTSFORENGLISH.csv';

fid = fopen(fname);

wordcounts = textscan(fid,'%s%s','Delimiter', ',');

% wordcounts = 

%     {97565x1 cell}    {97565x1 cell}
%   first cell is words  second is raw count of occurences



fclose(fid);

% convert the counts to numbers
rawcounts = str2num(char(wordcounts{2}));


% letters

letters = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N'...
    'O' 'P' 'Q' 'R' 'S' 'T' 'U' 'V' 'W' 'X' 'Y' 'Z'};


% not clear whether we care about doubled letters ('tt') or if pairs are
% ordered ('ab' vs 'ba').  this could matter for the independence analysis.
% for completeness we can count everything and then unpack

bigramcountmatrix = nan(26,26);


% as always, I'm guessing this is not very efficient

% let's do the matrix as rows

% pick a row
for i=1:26
    %     then for all letters
    for j=1:26
        % %         choose bigram
        letpair = strcat(letters(i), letters(j));
        fprintf('current bigram: %s\n',char(letpair));
        % %         count up the bigrams of up to n words
%         strmatch only finds strings that begin with the string
%         bigramcount=strfind(char(letpair),wordcounts{1}(1:nwords));
%         this works better but returns a cell
        bigramcount=strfind(wordcounts{1}(1:nwords),char(letpair));

%         %       the indices are given by 
%         find(cellfun(@sum,bigramcount))

%         you can get the list of words with the letterpair by
%         wordcounts{1}(find(cellfun(@sum,bigramcount)))

%         you can get the frequencies of the words with that letterpair by
        rawcounts(find(cellfun(@sum,bigramcount)));
    
%       there is an error to fix which is what happens when there are
%       mulitple instances of a pair within a word (like haha or
%       repetetive or teeter)

        


        %but we want to sum across all the occurences
        
        bigramcountmatrix(i,j) = sum(rawcounts(find(cellfun(@sum,bigramcount))));
 
    end

end


% visualize raw frequency matrix
figure('Name',['bigram counts for top ' num2str(nwords) 'English Words'],...
    'Color', [1 1 1],'Position',get(0,'ScreenSize'));

imagesc(bigramcountmatrix);
colorbar;
set(gca,'XTick',[1:26],'XTickLabel',letters,'YTick',[1:26],'YTickLabel',letters);


% now this doesn't tell us very much at least as a predictor.  for example
% if you want to use the fact that q is almost always followed by u as a
% measure, then you want to see the probability given q or each other
% letter.   so this just means normalizing the rows

% stupid way to do this
% for each row, divide it by its sum
for i=1:size(bigramcountmatrix,1)
    bigramcondprobmatrix(i,:)=bigramcountmatrix(i,:)/sum(bigramcountmatrix(i,:));
end



% visualize conditional probability  matrix
figure('Name',['conditional prob for top ' num2str(nwords) 'English Words'],...
    'Color', [1 1 1],'Position',get(0,'ScreenSize'));

imagesc(100*bigramcondprobmatrix);
h=colorbar;
caxis([0 10]);
set(gca,'XTick',[1:26],'XTickLabel',letters,'YTick',[1:26],'YTickLabel',letters);





% do probability preceded by
% stupid way to do this
% for each col, divide it by its sum
for i=1:size(bigramcountmatrix,1)
    bigramrevcondprobmatrix(:,i)=bigramcountmatrix(:,i)/sum(bigramcountmatrix(:,i));
end


% visualize conditional probability  matrix
figure('Name',['rev conditional prob for top ' num2str(nwords) 'English Words'],...
    'Color', [1 1 1],'Position',get(0,'ScreenSize'));

imagesc(100*bigramrevcondprobmatrix);
h=colorbar;
caxis([0 10]);
set(gca,'XTick',[1:26],'XTickLabel',letters,'YTick',[1:26],'YTickLabel',letters);







% you could go one step deeper and ask which letters are similar in terms
% of the letters they proceed
figure;
plotmatrix(bigramcondprobmatrix');

bigramcontextsim = corr(bigramcondprobmatrix','rows','pairwise');

figure('Name',['similar contexts for brigram cond prob ' num2str(nwords) 'English Words'],...
    'Color', [1 1 1],'Position',get(0,'ScreenSize'));

imagesc(bigramcontextsim);
h=colorbar;
% caxis([0 10]);
set(gca,'XTick',[1:26],'XTickLabel',letters,'YTick',[1:26],'YTickLabel',letters);

