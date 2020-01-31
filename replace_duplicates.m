function newLine = replace_duplicates(line,T)
    newLine = line;
    duplicatesIndices = logical([0 ~diff(line)]);
    yetIn = line(~duplicatesIndices);
    addable = setdiff(randperm(T),yetIn,'stable');
    newLine(duplicatesIndices) = addable(1:sum(duplicatesIndices));
end
