function X = normalizeM(x)
    ave = mean(x,2);
    x = x - ave*ones(1,size(x,2));
    s = mean(abs(x),2);%std(x,0,2);
    X = x ./ (s*ones(1,size(x,2)));
end