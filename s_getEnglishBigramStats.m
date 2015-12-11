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
nwords = 50000;


% load the data

fname = 'WORDCOUNTSFORENGLISH.csv';

fid = fopen(fname);

wordcounts = textscan(fid,'%s%s','Delimiter', ',');

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
        
        bigramcount=strmatch(char(letpair),wordcounts{1}(1:nwords));
%         
%         have to multiply instances by count of the word
        


        %             then sum all found instances
        
        bigramcountmatrix(i,j) = sum(cellfun(@sum,bigramcount));
 
    end

end


% visualize matrix
figure('Name',['bigram counts for top ' num2str(nwords) 'English Words'],...
    'Color', [1 1 1],'Position',get(0,'ScreenSize'));

imagesc(bigramcountmatrix);
colorbar;
set(gca,'XTick',[1:26],'XTickLabel',letters,'YTick',[1:26],'YTickLabel',letters);




