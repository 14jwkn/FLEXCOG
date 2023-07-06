Adapts code from: 
https://github.com/trendscenter/gift

function [bestx, yfit] = fit_L_tocurve_area(x,y)
%[bestx, F] = fit_L_tocurve_area(x,y, plotit)
%
% if nargin < 3
%     plotit = 1;
% end
x = x(:);
y = y(:);
%%
% P = [xbp m1 m2 b]
options = optimset('TolFun',1e-14, 'TolX', 1e-14, 'MaxFunEvals', 100000, 'MaxIter', 10000);
P0(1) = x(5);
P0(2) = (y(4)-y(2))/(x(4)-x(2));
P0(3) = (y(end)-y(end-1))/(x(end)-x(end-1));
P0(4) = y(5);

LB = [x(2)    -Inf -Inf  min(y)];
UB = [x(end-1) 0     0   max(y)];

[PF,RESNORM,RESIDUAL,EXITFLAG,OUTPUT] = icatb_lsqcurvefit(@Lcurve,P0,x,y, LB, UB, options);

bestx = ceil(PF(1));

yfit = Lcurve(PF, x);

function Z = Lcurve(P, XDATA)
% P = [xbp m1 m2 b]

XDATA = XDATA-P(1);

% curve 1
Y1 = P(2)*XDATA + P(4);
% curve 2
Y2 = P(3)*XDATA + P(4);

Z = zeros(size(Y1));
Z(XDATA < 0) = Y1(XDATA < 0);
Z(XDATA >= 0) = Y2(XDATA >= 0);


