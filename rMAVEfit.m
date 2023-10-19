function [yh,yp,B0,hm0,cv0,h0] = rMAVEfit(X,y,m0,Xp,h0,flag)
% Automatic estimation of the rMAVE model
% Reqire function "rMAVE.m", "ksLLadap.m", "cvadap.m"
% Tested on Matlab R2010a
% Input:
%     X --- n x p matrix for expanaltory variables
%     y --- n x 1 vector for response variable
%     Xp --- np x p matrix for out-of-sample predictions
%     m0 --- dimension of the effective DRS
% Output:
%     B0 --- p x m0 matrix with each column be the estimated piecewise 
%            single-index
%     yh --- n x 1 vector for fitted response variable
%     yp --- np x 1 vector for predictions at Xp
%     h0 --- 1 x 1 scalar for the bandwidth used in kernels
% Reference: Y. Xia, H. Tong, W.K. Li and L.X. Zhu,
% "Adaptive estimation of the effictive dimension space",(2002)
% ---------------------------------------------------------------------- %
if nargin < 4
    Xp = [];
end
if nargin < 6
    flag = 1;
end
cv0 = 1e3 * sum(abs(y));
n = size(X,1);

if nargin < 5 || isempty(h0)
    hrot = n^(-1 / (m0 + 4));
    htot = (1 : .4 : 1.8) * hrot;
else
    htot = h0;
end

for h = htot
    B = rMAVE(X,y,h,m0);
    B = real(B);
    [hm,fm,cv] = cvadap(X * B,y);
    if cv < cv0
        cv0 = cv;
        h0 = h;
        hm0 = hm;
        f0 = fm;
        B0 = B;
    end
end
mf = mean(f0);

yh = [];
yp = [];
if flag
    XB = X * B0;
    yh = ksLLadap(XB,y,XB,hm0,mf);
    if ~isempty(Xp)
        yp = ksLLadap(XB,y, Xp * B0, hm0,mf);
    end
end


