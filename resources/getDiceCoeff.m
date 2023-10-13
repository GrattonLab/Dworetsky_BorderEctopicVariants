function DC = getDiceCoeff(A,B)
% 
% call:
% 
%        DC = getDiceCoeff(A,B)
%        
% Compute the Dice Coefficient, a measure of similarity between two binary (0,1) vector-sets A, B.
% E.g. to compute the DC between two network community partitions, first assign each
% link to the corresponding community (e.g. using 'getCommunityMatrix.m'), then binarize the cor-
% responding matrix, extract the sub-diagonal elements in form of a vector. These vectors are the
% required inputs, A and B.
% 
% Between others, alternative measures of similarity between sets are:
% 
% - Normalized Mutual Information ('getNMI.m') 
% - Jaccard Index ('getJaccard.m')
% 
% INPUT
% 
%        A  :  Binary vector of set A
%        B  :  Binary vector of set B
%        
% OUTPUT
% 
%        DC  :  Dice Coefficient
%        
%        
% R.G. Bettinardi, Computational Neuroscience group, UPF, Barcelona, Spain, (rug.bettinardi@gmail.com), April 2016
% ----------------------------------------------------------------------------------------------------------------

if nargin<2
    error('The function requires two input-vectors')
end

if (sum(A==0) + sum(A==1)) < length(A)
    error('Input A must be binary (0,1) !!!')
end
  
if (sum(B==0) + sum(B==1)) < length(B)
    error('Input B must be binary (0,1) !!!')
end
       
DC = 2*(sum(A.*B))/sum(A + B);