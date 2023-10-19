function [B, WB, a, ab, f] = rOPGadapFit(xx, y, h0, nd, B0, wt0)
%Searching the effective dimension reduction subspace of model
%  			y = g(B^TX) + e
%Useage: directions(x, y, h, d, g)
%Input: 
%     x --- expanaltory variables
%     y --- response
%     h --- bandwidth
%     nd --- working dimension of the space
%     B0 --- user supplied initial value, default to be eye(p)
%     wt0 --- user supplied weights for weighted average OPG
%Output:
%     B --- e.d.r directions
%     WB --- weights for B
%     a --- fitted y
%     ab ---  fitted pointwise gradient
%     f --- nonparametric density of X
%Reference: Y. Xia, H. Tong, W.K. Li and L.X. Zhu, 
% "Adaptive estimation of the effictive dimension space",(2002)  
% ------------------------------------------------------------
[n,p] = size(xx);
mm = mean(xx);
xx = xx - repmat(mm,n,1);
ss = pinv(xx'*xx/n)^0.5;
xx = xx*ss;

if  nargin < 3 || isempty(h0)
    h0 = 1.06 * n^(-1 / (nd + 4))  / (1 + log(nd)); %%%%%%%%%%%%%%%
end
if nargin < 4 || isempty(nd)
    nd = p;
end
if nargin < 5 || isempty(B0)
    B0 = eye(p);
    B = B0;
    flagB = 1;
else
    B = pinv(ss) * B0;
    flagB = 0;
end

if nargin < 6 || isempty(wt0)
    wt0 = ones(n,1);
end
wt0 = wt0 / sum(wt0);

% Find the proper initial h
M_Sdist = zeros(n,n);
XB = xx * B;
Sdist = pdist(XB);
count = 0;
for i = 1:(n - 1)
    M_Sdist((i+1):end,i) = Sdist((1+count):(n-i+count));
    count = count + n - i;
end
M_Sdist = M_Sdist + M_Sdist';
npk = max(p + 1,round(sqrt(n / nd)));
hdist = zeros(n,size(B,2));
for j = 1 : n
    Xdist = M_Sdist(:,j);
    [~,Ij] = sort(Xdist);
    Ij_t = Ij(1 : npk,:);
    XXdist = ones(npk,1) * XB(j,:) - XB(Ij_t,:);
    hdist(j,:) = max(abs(XXdist));
end
hdist = median(hdist);
% ------------------------- %

onen = ones(n,1);
Btemp = zeros(p,nd);
f = ones(n,1);
WB = ones(size(B,2),1);
cf2 = 1.2;
for iter = 1 : 10 + 10 * flagB;
    ab = ones(p,n);
    a = y;    
    XB = xx * B;
    
    hp = size(B,2);
    
    if iter <= 1 % && size(B,2) >= p / 2
        h = hdist;
    else
        h = max(h0 * std(XB) ./ vec(WB)', ...
            hdist * (B.^2) * n^(-(iter - 1) / (hp + 6) / 2));
    end
    
    for i = 1 : n;
        xiB = XB - repmat(XB(i,:),n,1);
        H = ones(n,1) * h;
        xd = xiB ./ H;
        kernel = .75^hp * prod((1 - xd.^2).*(abs(xd)<1) ./ H, 2);
        f(i) = mean(kernel);
    end
    maxf = max(f);
    f = f + max(f) / n;
    f2 = (mean(f) ./ f).^cf2;
    for i = 1 : n;
        xi = xx - repmat(xx(i,:),n,1);
        xiB = XB - repmat(XB(i,:),n,1);
        
        h2 = max([h * min(f2(i),(sum(WB) + 1)^cf2); h]); 
        H2 = ones(n,1) * h2;
        xd = xiB ./ H2;
        kernel = prod(exp(-xd.^2 / 2) ./ H2, 2) / sqrt(2 * pi)^hp + min(maxf,1) / (n * n);
        
        onexi = [onen xi];
        xk = onexi.*repmat(kernel, 1, p+1);
        XWY = xk' * y;
        abi = pinv(xk'*onexi)*XWY;
        ab(:,i) = abi(2:p+1);
        a(i) =  abi(1);        
    end;
    normab = sum(ab.^2);
    qab1 = quantile(normab,.25);
    qab3 = quantile(normab,.75);
    ubnd = qab3 + 1.5 * (qab3 - qab1);
    wab = (normab < ubnd);
    normab2 = normab .* wab + ubnd * (1 - wab);
    ab = ab .* (ones(p,1) * (normab2 ./ normab));
    
    qf1 = quantile(f,.25);
    qf3 = quantile(f,.75);
    fqw = (f > qf1 - 1.5 * (qf3 - qf1));
    
    ay = (a - y).^2;
    qay1 = quantile(ay,.25);
    qay3 = quantile(ay,.75);
    aerr = (ay < 1.5 * (qay3 - qay1) + qay3); %
    fw = sqrt(fqw .* wt0 .* aerr);
    Mab = ab * diag(fw) *ab' / (sum(fw) + eps);
    ssMab = sum(sum(Mab));
    if isnan(ssMab) || isinf(ssMab)
        B = Btemp;
        Ds = zeros(p,1);
    else
        [B0 D] = eig(Mab);
        [D, I] = sort(abs(diag(D)),'descend');        
        Itemp = I(1 : nd);
        B = B0(:,Itemp);
        Dg = abs(D(1 : nd));
        Ds = cumsum(Dg / sum(abs(D)));
    end
    if norm(B - Btemp) < 2 * sqrt(p) * 1e-5
        break
    else
        Btemp = B;
    end
    flagDs = (Ds > .99);
    [~,iDs] = max(flagDs);
    if iDs < nd
        WB = ones(nd,1) * eps;
        WB(1:max(iDs,min(2,nd))) = 1;
    else
        WB = ones(nd,1);
    end
end;

B = real(ss*B);
for i = 1:size(B,2);
    B(:,i) = B(:,i)/sqrt(B(:,i)'*B(:,i));
end
ab = real(ab'*ss');