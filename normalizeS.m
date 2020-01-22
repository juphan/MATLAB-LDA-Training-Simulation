function X = normalizeS(x)
    ave = mean(x,2);
    x = x - ave*ones(1,size(x,2));
    s = std(x,0,2);
    X = x ./ (s*ones(1,size(x,2)));
end