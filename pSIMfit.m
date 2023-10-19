function [pBfit0,yh,yp,gI0,gIp,cv,BIC,Bw0,mf0,fp,Btot0,h0]...
    = pSIMfit(X,y,m0,Xp,hrot)
% Automatic estimation of the piecewise single index model
% Reqire function "pSIM.m","rMAVEfit", "ksLLadap.m", "cvadap.m"
% Tested on Matlab R2010a
% Input:
%     X --- n x p matrix for expanaltory variables
%     y --- n x 1 vector for response variable
%     m0 --- number of effective pieces
%     Xp --- np x p matrix for out-of-sample predictions
%     hrot --- bandwiths to be selected
%     m1 --- the effective dimension
% Output:
%     pB --- p x m0 matrix with each column be the estimated piecewise 
%            single-index
%     yh --- n x 1 vector for fitted response variable
%     yp --- np x 1 vector for predictions at Xp
%     gI0 --- n x 1 vector for estimated group index
%     h0 --- m0 x 1 vector of bandwidths used in kernel smoothing for 
%             each region
%     t0 --- 1 x 1 scalar for the number of neighors used for local
%            gradients estimation
% Remark: "h0" is for calculation of BIC.
% Reference: Wang, T. and Xia, Y.,
%  "A piecewise single index model for dimension reduction",(2012)
% ---------------------------------------------------------------------- %
if nargin < 4
    Xp = [];
end
[n,p] = size(X);

cv0 = 1e3 * var(y);
h = zeros(m0,1);
mf = zeros(m0,1);
f = zeros(n,1);

if nargin < 5 || isempty(hrot)
    hrot = n^(-1 / (m0 + 4)) / (1 + log(m0));
%     hrot = [1.2,1.5,1.8] * n^(-1 / (m0 + 4)) / (1 + log(m0));
end

% mm = mean(X);
% Xm = X - repmat(mm,n,1);
% iss = (Xm' * Xm / n)^0.5; % pinv(ss)

% mm = mean(X);
% X = X - repmat(mm,n,1);
% ss = pinv(X'*X/n)^0.5;
% X = X*ss;

Htot = hrot; %%%%%%%%%%%%%%%%%%%%%%%%
nHtot = length(Htot);
cv = zeros(nHtot,1);
ict = 1;

pw = min(m0,p);
sqrn = sqrt(n);
for hOPG = Htot
    BwI = [];
    for iter = 1 : 1 
        if m0 == 1 % Single-index
            pB = rMAVE(X,y,.5,1);
            Bw = pB;
            Wb = 1;
            pBfit = pB;
            gI = ones(n,1);
            IBtot = ones(n,1);
            fBw = ones(n,1);
            Btot = ones(n,p);
        else
            [Bw, Wb, ~, Btot, fBw,IBtot] = rOPGadapH(X,y,hOPG,pw,BwI);
%             Bw = rMAVE(X,y,.5,size(Bw,2));%%%%%%%%%%%%%%%%%
            if isempty(Xp)
                IBtot = ones(n,1);
            end
            [pB,gI,Bpnt] = pSIMT(X,y,m0,Bw,Wb,Btot,IBtot);
            pBfit = pB;
            
            ng = zeros(m0,1);
            
            for mp = 1 : m0
                Imp = find(gI == mp);
                ng(mp) = length(Imp);
                if ng(mp) < m0
                    ng(mp) = ng(mp) + eps;
                    continue
                end
                
%                 if p >= sqrn%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     pBtemp = pB(:,mp) .* diag(iss);
%                     pBtemp = pBtemp / norm(pBtemp);
%                     pB(:,mp) = pB(:,mp) .* (abs(pBtemp) > 1 / sqrn);
%                     pB(:,mp) = pB(:,mp) / norm(pB(:,mp));
%                 end
            end
        end
%         pBpB = pB*diag(ng)*pB'/sum(ng) + Bw*diag(Wb)*Bw'/sum(Wb);
%         [V,D] = eig(pBpB);
%         [D,I] = sort(abs(diag(D)),'descend');
%         D = D / sum(D);
%         nd = min(max(sum(D >= .05),2),m0);
%         BwI = V(:,I(1:nd));
    end
    % + CV prediction at each point based on the predicted group    
    for m = 1 : m0
        Im = find((gI == m).*IBtot);
        nm = length(Im);
        if nm < 1
            continue
        end
        
        XBm = X(Im,:) * pB(:,m); %%%%%%%%%%%%%%%%%%
        ym = y(Im,:);
        [hm,fm,cvm] = cvadapW(XBm,ym);
        mf(m,:) = mean(fm);
        f(Im,:) = fm;
        h(m,:) = hm;        
        cv(ict,:) = cv(ict,:) + cvm;
    end
    % ----------------------------------------------- %
    if cv(ict,:) < cv0
        cv0 = cv(ict,:);
        pB0 = pB;
        pBfit0 = pBfit;
        gI0 = gI;
        h0 = h;
        hOPG0 = hOPG;
        fB0 = fBw;
        mf0 = mf;
        Bw0 = Bw;
        Wb0 = Wb;
        Btot0 = Btot;
        f0 = f;
    end
    ict = ict + 1;
end
% In-sample fitting and Out-of-sample prediction
yh = zeros(n,1);
np = size(Xp,1);
yp = zeros(np,1);
fp = zeros(np,1);
ng = ones(m0,1);
XB = X * Bw0;
stdXB = min(std(XB),iqr(XB)/1.34);
pBw0 = size(Bw0,2);
h2 = n^(-1/(4 + pBw0)) * stdXB; %hOPG0 * (stdXB ./ vec(Wb0)');
H2 = ones(n,1) * h2;

IRB = 1:n;%find(IBtot);
if np > 0 && m0 > 1
    %%%%%%%%%%%%%%%%%%%%    
    XRB = X(IRB,:);
    gIRB = gI0(IRB);
    %%%%%%%%%%%%%%%%%%%%
%     gIp = classify(Xp * Bw0,XRB * Bw0,gIRB,'quadratic');
%     Bdry = zeros(n,1);
%     knnWt = zeros(n,m0);
    logn = fix(log(n)) + 1;
    maxfB0 = min(quantile(fB0,.99),(2*pi)^(-pBw0/2));    
    [gIp, Bdry, knnWt] = knn(XRB,gIRB,Xp,logn, Bw0 * diag(sqrt(Wb0))); 
    fBdry = ones(np,1);
    for i = 1 : np
        XBi = Xp(i,:) * Bw0;
        xd = (XB - ones(n,1) * XBi) ./ H2;
        xd2 = xd.^2;
        fBdry(i) = mean(prod(exp(-xd2 / 2) ./ H2,2)) / sqrt(2 * pi)^pBw0;
    end
    % Sparse points
    flagfB = (fBdry < maxfB0 * exp(-2));
    ISps = find(flagfB);
    if ~isempty(ISps)
        [gIp(ISps,:), Bdry(ISps,:), knnWt(ISps,:)] = ...
            knn(XRB,gIRB,Xp(ISps,:),logn * 2, Bw0 * diag(sqrt(Wb0)));
    end
elseif np > 0
    gIp = ones(np,1);
    Bdry = zeros(np,1);
    knnWt = ones(np,1);
else
    gIp = [];
    yp = [];
    Bdry = [];
    knnWt = [];
end

% for m = 1 : m0    
%     Ipm = find(gIp == m);
%     Im = find(gI0 == m);
%     XBRB = X(Im,:) * Bw0;
%     mXB = repmat(max(XBRB),length(Ipm),1);
%     Xpm =  Xp(Ipm,:) * Bw0;
%     Imm = find(sum(Xpm > mXB,2));
%     Bdry(Ipm(Imm),:) = 1;
%     knnWt(Ipm(Imm),:) = 1;
% end
sum(Bdry)
%%%%%%%%%%%%%%

stdXpB = ones(m0,1);
% ASE = 0;
% ASEgp = ones(m0,1) * std(y) * 1e3;
Cntr = zeros(m0,size(Bw0,2));
Ibd = zeros(n,1);
for m = 1 : m0
    Im = find(gI0 == m);
    
    nm = length(Im);
    if nm < log(n)
        ng(m) = nm + eps;
        continue
    end
    
    ym = y(Im,:);
    pBm = pB0(:,m);
    XBm = X(Im,:) * pBm;
    Ibd(Im,:) = max((XBm > quantile(XBm,.99)),(XBm < quantile(XBm,.01)));
    f0m = f0(Im,:);
    stdXpB(m) = min(std(XBm),iqr(XBm)/1.34);
    ng(m) = nm;
    Cntr(m,:) = median(XB(Im,:) * diag(Wb0));
    
    yh(Im,:) = ksLLadapW(XBm,ym,XBm,h0(m),mf0(m));     
end
% % Outliers
% ay = (yh - y).^2;
% may = mean(ay);
% ubd1 = chi2inv((n - m0)/n,fix(may)+1);
% ubd2 = quantile(ay,(n - m0)/n);
% fay = (ay > max(ubd1,ubd2));
% ifay = find(fay);
% nay = length(ifay);
% if nay > 0
%     Xay = X(ifay,:);
%     yay = y(ifay,:);
%     
%     yayh = zeros(nay,m0);
%     for m = 1 : m0
%         Im = find(gI0 == m);
%         nm = length(Im);
%         if nm < log(n)
%             ng(m) = nm + eps;
%             continue
%         end
%         
%         ym = y(Im,:);
%         pBm = pB0(:,m);
%         XBm = X(Im,:) * pBm;
%         yayh(:,m) = ksLLadapW(XBm,ym,Xay*pBm,h0(m),mf0(m));
%     end
%     [~,gI0(ifay,:)] = min(abs(yayh - repmat(yay,1,m0)),[],2);
%     yh(ifay,:) = reshape(diag(yayh(1:nay,gI0(ifay,:))),nay,1);
% end
ASEtot = (y - yh).^2;
ASE = sum(ASEtot);

% Predictions
% NIRB = find(~IBtot);
% y(NIRB,:) = yh(NIRB,:);
if ~isempty(Xp)
    for m = 1 : m0
        Im = find((gI0 == m) .* max(IBtot,Ibd));
        nm = length(Im);
        if nm < log(n)
            ng(m) = nm + eps;
            continue
        end
        
        ym = y(Im,:);
        pBm = pB0(:,m);
        XBm = X(Im,:) * pBm;      
        
        Ipm = find(gIp == m);
        Ipm = setdiff(Ipm,find(Bdry)); % The inner points of Xp
        npm = length(Ipm);
        
        if npm > 0
            Xpm = Xp(Ipm,:) * pBm;
            
            ff = zeros(npm,1);
%             h = max(nm^(-0.2) * stdXpB(m) / 2, h0(m));
%             H = ones(nm,1) * h;
%             
%             for i = 1 : npm
%                 xd = (XBm - ones(nm,1) * Xpm(i,:)) ./ H;
%                 xd2 = xd.^2;
%                 ff(i) = mean(.75 * prod((1 - xd2).*(xd2 < 1) ./ H,2));
%             end
%             ff = ff + mf0(m) / n;
            fp(Ipm) = ff;
            yp(Ipm,:) = ksLLadapW(XBm,ym,Xpm,h0(m),mf0(m));
        end
    end
end
hsmooth = h0;%max([h0,ng.^(-0.2) .* stdXpB],[],2);
BIC = log(ASE / n) + log(n) * sum(1 ./ ng ./ hsmooth .* stdXpB);
%%%%%%%%%%%%%%%%
IBdry = find(Bdry);
nBdry = length(IBdry);
sASEgp = ones(m0,1);% sqrt(ASEgp + min(ASEgp(ASEgp > 0)));
if nBdry > 0 % The boundary points in Xp
    WtBy = knnWt(IBdry,:) .* (ones(nBdry,1) * (1 ./ sASEgp'));    
    ypBdry = zeros(nBdry,m0);
    for m = 1 : m0
        Im = find((gI0 == m) .* max(IBtot,Ibd));
        if isempty(Im)
            continue
        end
        nm = length(Im);
        ym = y(Im,:);
        pBm = pB0(:,m);
        XBm = X(Im,:) * pBm;
        ng(m) = nm;
        
        Ipm = IBdry;
        npm = length(Ipm);        
        Xpm = Xp(Ipm,:) * pBm;
        gIpm = gIp(Ipm,:);
        
%         ff = zeros(npm,1);
%         h = max(nm^(-0.2) * stdXpB(m) / 2, h0(m));
%         H = ones(nm,1) * h;
%         
%         for i = 1 : npm
%             xd = (XBm - ones(nm,1) * Xpm(i,:)) ./ H;
%             xd2 = abs(xd);
%             ff(i) = sum((xd2 < 1));%mean(.75 * prod((1 - xd2).*(xd2 < 1) ./ H,2));
%         end
%         fp(Ipm) = ff;
%         Ipmt = find((ff > 0) + (gIpm == m));
%         if ~isempty(Ipmt)
%             ypBdry(Ipmt,m) = ksLLadapW(XBm,ym,Xpm(Ipmt,:),h0(m),mf0(m));
%         end
%         NIpmt = setdiff(1:npm,Ipmt); 
%         if ~isempty(NIpmt)
%             WtBy(NIpmt,m) = 0;
%         end
        ypBdry(:,m) = ksLLadapW(XBm,ym,Xpm,h0(m),mf0(m));
    end    
    WtBy = WtBy ./ repmat(sum(WtBy,2),1,m0);
    yp(IBdry,:) = sum(ypBdry .* WtBy,2);
end

% pBfit0 = real(ss*pBfit0);
% for i = 1:m0
%     pBfit0(:,i) = pBfit0(:,i)/sqrt(pBfit0(:,i)'*pBfit0(:,i));
% end

% nISps = length(ISps);
% % ISps, Xp(ISps,:) * Bw0, Xp(ISps,:) * pB0,
% if nISps > 0 % The sparse points in Xp
%     ypb = zeros(nISps,m0);      
%     DistC = pdist(Cntr);
%     minC = min(DistC);
%     for i = 1 : nISps
%         mi = gIp(ISps(i));        
%         XBi = Xp(ISps(i),:) * Bw0 * diag(Wb0);
%         Disti = sqrt(sum((ones(m0,1) * XBi - Cntr).^2,2));
%         [~,IDisti] = sort(Disti);
%         min1 = Disti(IDisti(1));
%         min2 = Disti(IDisti(2));
%         if (min1 + minC / 2) < min2
%             ypb(i,:) = yp(ISps(i),:) * mean(1 ./ sASEgp);
%             continue
%         else
%             ypb(i,mi) = yp(ISps(i),:) / sASEgp(mi);
%         end
%         for m = setdiff(1 : m0, mi)
%             Im = find(gI0 == m);
%             if isempty(Im)
%                 continue
%             end
%             nm = length(Im);
%             ym = y(Im,:);
%             pBm = pB0(:,m);
%             XBm = X(Im,:) * pBm;
%             Xpi = Xp(ISps(i),:) * pBm;
%             
%             h = max(nm^(-0.2) * stdXpB(m) / 2, h0(m));
%             H = ones(nm,1) * h;
%             xd = (XBm - ones(nm,1) * Xpi) ./ H;
%             xd2 = xd.^2;
%             ffi = mean(.75 * prod((1 - xd2).*(xd2 < 1) ./ H,2)); 
%             flagffi = (ffi > min(1 / (nm * h), 2 / nm^.8));            
%                         
%             if flagffi
%                 ypb(i,m) = ksLLadap(XBm,ym,Xpi,h0(m),mf0(m)) / sASEgp(m);
%             else
%                 ypb(i,m) = yp(ISps(i),:) / sASEgp(m);
%             end
%         end        
%     end
%     yp(ISps,:) = sum(ypb,2) / sum(1 ./ sASEgp);
% end
% ---------------------------------------------------------------------- %

% May exclude the training points with extremely low densities
