function [yh,yp,B0,Dg0,Btot0,hm0,h0] = rOPGfit(X,y,m0,Xp,flag)
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
if nargin < 5
    flag = 1;
end
cv0 = 1e3 * sum(abs(y));
n = size(X,1);
hrot = 1.06 * n^(-1 / (m0 + 4))  / (1 + log(m0));

for h = hrot % (.8 : .2 : 1.2) * hrot
    [B, Dg, a, Btot] = rOPGadap(X,y,h,m0);
%     Dxb = sqrt(Dg) / max(sqrt(Dg)) + 1 / n^2;
    [hm,fm,cv] = cvadap(X * B,y); % ,[],[],Dxb);
    if cv < cv0
        cv0 = cv;
        h0 = h;
        hm0 = hm;
        f0 = fm;
        B0 = B;
        Dg0 = Dg;
        Btot0 = Btot;
    end
end
mf = mean(f0);

yh = [];
yp = [];
if flag
    XB = X * B0;
    yh = ksLLadap(XB,y,[],hm0,mf);
    if ~isempty(Xp)
        yp = ksLLadap(XB,y, Xp * B0, hm0,mf);
    end
end

