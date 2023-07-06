Adapts code from: 
https://github.com/trendscenter/gift

%Define command line arguments.
function [results] = elbow_k(X,klist,clusts,cents)

%Collect WCS and CVI.
SSE = zeros(1, length(klist));
R = zeros(1,length(klist));

%Evaluate clustering performance.
for i=1:size(klist,2)
    
    %Set up values.
    idx = clusts{i};
    C = cents{i};
    D = distfun(X,C,'cityblock',1,1,1);
    n = size(idx,1);
    d = D((idx-1)*n + (1:n)');
    sumD = accumarray(idx,d,[klist(i),1]);
    
    %Elbow method.
    [~, R(i)] = cluster_goodness(D,idx);
    tmpSumD = sumD;
    SSE(i) = sum(tmpSumD(:));
end

%Find optimal k.
[w_bestx,w_yfit] = fit_L_tocurve_area(klist,SSE);
[c_bestx,c_yfit] = fit_L_tocurve_area(klist,R);
results.wcs = SSE;
results.cvi = R;
results.w_y = w_yfit;
results.w_k = w_bestx;
results.c_y = c_yfit;
results.c_k = c_bestx;
end
