function [CPM] = getCondPMatrix(data)

% function [CPM] = getCondPMatrix(data)
% makes a conditional probability matrix CPM from the the matrix data
% data is n subjects x 26 letters where each entry is a color label
% cpm is a conditional probability matrix where each row can be interpreted
% as given that a letter was assigned a color, what is the probability of
% all the other possible letter color assignments



% make a large single figure with all the letters vs all the letters so
%make indices so we know where to place each probability distribution we
%solve for
xindvector = ((0:26)*12)+1;
% matrix to hold our solutions
CPM = zeros(12*26);
% for each letter
for lA=1:26
    %     make a variable to store the probability distributions
    
    %     pick another letter (only need upper triangle
    for lB=1:26
        %     for each color category
        for cA=1:12
            %     get a histogram of the values of lB
            %     get index to entries where letter A has color a
            lAcA = find(data(:,lA)==cA-1);
            %          then get all of the values of B for those subjects
            temp = hist(data(lAcA,lB),0:11);
            %             turn this into a probability distribution?
            temp = temp/sum(temp);
%             inset into our matrix
            CPM(xindvector(lA)+cA-1,xindvector(lB):xindvector(lB)+11)=temp;
        end

    end
end






end