function r = cell2sym(c)
c = c(:);
r = sym(zeros(numel(c),1));
for I=1:length(c)
    r(I) = c{I};
end
