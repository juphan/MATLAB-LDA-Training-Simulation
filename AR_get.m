function [Feat] = AR_get(x,wl,wi,order)
% Make sure your input for x is cut_data'

datasize = size(x,2);
Nsignals = size(x,1);
numwin = floor((datasize-wl)/wi);

Feat = [];

for f=1:order+1
    Feat(f).data = zeros(Nsignals, numwin);
    Feat(f).name = sprintf('AR(n-%d)', f);
end
Feat(f).name = 'AR(err)';

strt = 1+wi;
endn = wl+wi;

for w = 1:numwin
    n = 0;
    for i=1:Nsignals      %vert
        n = n+1;
        cx = x(i, strt:endn);
        for dt = 1:order
            cx = vertcat(cx,x(i, strt-dt:endn-dt));
        end
        cX = normalizeS(cx);
        B = regress(cX(1,:)',cX(2:end,:)');
        y = cX(2:end,:)'*B;
        error = cX(1,:)'-y;
        err = (sum(error.^2)/wl)^0.5;
        for f = 1:order
            Feat(f).data(i,w) = B(f);
        end
        Feat(f+1).data(i,w) = err;
    end
    strt = strt+wi;
    endn = endn+wi;
end

end