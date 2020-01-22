function prediction = LDA_Predict(features, Wg, Cg)

tmp = features'*Wg + Cg;

max = tmp(1);
gtemp = 1;

for i=2:5
    if (tmp(i)>max)
        max = tmp(i);
        gtemp = i;
    end
end

prediction = gtemp;

end