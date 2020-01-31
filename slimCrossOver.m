function [firstKid,secondKid] = slimCrossOver(firstParent,secondParent, firstKid,secondKid)
% For two parents A and B (firstParent and secondParent), this function compute
% the part of the DNA of parents A and B that are in the sets A\B and in
% B\A. Then a subset A' of A\B and a subset B' of B\A are swapped (note
% that these sets have the same cardinal, i.e. |A'|==|B'|). It means that
% new A:= (A\A') union B'
% new B:= (B\B') union A'

% Find part of the DNA of parents A and B that are in A\B and in
% B\A.
[~,indFirstParent,indSecondParent] = setxor(firstParent,secondParent,'sorted');

% random shuffle
lenIndParent = length(indFirstParent);
if lenIndParent ~= 0
    indFirstParent = indFirstParent(randperm(lenIndParent));
    indSecondParent = indSecondParent(randperm(lenIndParent));
    
    numberOfSwaps = unidrnd(lenIndParent);
    % keep only some swaps
    indFirstParent = indFirstParent(1:numberOfSwaps);
    indSecondParent=indSecondParent(1:numberOfSwaps);
    
    firstKid(indFirstParent) = secondParent(indSecondParent);
    secondKid(indSecondParent)=firstParent(indFirstParent);
end
end