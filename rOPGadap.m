function [B, WB, a, ab, fp, Ifw] = rOPGadap(xx, y, nd, B0, wt0)
%Searching the effective dimension reduction subspace of model
%  			y = g(B^TX) + e
%Useage: directions(x, y, h, d, g)
%Input:
%     xx --- expanaltory variables
%     y --- response
%     nd --- working dimension of the space
%     B0 --- user supplied initial value, default to be eye(p)
%     wt0 --- user supplied weights for weighted average OPG
%Output:
%     B --- e.d.r directions
%     WB --- weights for B
%     a --- fitted y
%     ab ---  fitted pointwise gradient
%     fp --- nonparametric density of X
%     Ifw --- index of outliers
%Reference: Y. Xia, H. Tong, W.K. Li and L.X. Zhu,
% "Adaptive estimation of the effictive dimension space",(2002)
% ------------------------------------------------------------
[n,p] = size(xx);
mm = mean(xx);
xx = xx - repmat(mm,n,1);
ss = real(pinv(xx'*xx/n)^0.5);
xx = xx*ss;

if nargin < 3 || isempty(nd)
    nd = p;
end
if nargin < 4 || isempty(B0)
    B0 = rMAVE(xx,y,.5,nd);
    B = B0;
else
    B = pinv(ss) * B0;
end

if nargin < 5 || isempty(wt0)
    wt0 = ones(n,1);
end
wt0 = wt0 / sum(wt0);

Btemp = zeros(p,p);
Btemp2 = Btemp;
WB = ones(size(B,2),1);
Ifw = ones(n,1);
fp = ones(n,1);
ITER = 5 * nd + 5;
for iter = 1 : ITER
    Idit = find(Ifw);
    nfw = length(Idit);
    onen = ones(n,1);
    ab = ones(p,n);
    a = ones(n,1);
    cviter = ones(n,1);
    xxit = xx(Idit,:);
    yit = y(Idit,:);
    XB = xx * B;
    XBWB = XB * diag(sqrt(WB));
    XBWBit = XBWB(Idit,:);
    sigXB = min(iqr(XB)/1.34,1);
    
    hp = size(B,2);
    HCV = n^(-1/(hp+4))*(2:.5:(log(n)+1)^1.5)/log(n)*log(p/2);
    Lhcv = length(HCV);
    Licv = zeros(Lhcv,1);
    CVtot = zeros(n,Lhcv);
    ABtot = zeros(p,n,Lhcv);
    Atot = zeros(n,Lhcv);
    FPtot = zeros(n,Lhcv);
    KERtot = zeros(n,n,Lhcv);
    ic = 1;
    for hcv = HCV
        h = hcv * sigXB;
        H = onen * h; 
        spihp = (2 * pi)^(hp/2);
        
        for i = 1 : n;
            if Ifw(i)
                ii = sum(Ifw(1:i));
                xi = xxit - repmat(xxit(ii,:),nfw,1);
                xiB = XBWBit - repmat(XBWBit(ii,:),nfw,1);
                yi = yit;
                Hi = H(1:nfw,:);
                xd = xiB ./ Hi;                
                Iditi = Idit;
            else
                ii = 1;
                xi = xxit - repmat(xx(i,:),nfw,1);
                xi = [zeros(1,p);xi];
                xiB = XBWBit - repmat(XBWB(i,:),nfw,1);
                xiB = [zeros(1,hp);xiB];
                yi = [y(i,:);yit];
                Hi = H(1:nfw+1,:);
                xd = xiB ./ Hi;
                Iditi = [i;Idit];
            end            
            
            xd2 = xd.^2;
            kernel = prod(exp(-xd2 / 2) ./ Hi, 2) / spihp + eps;
            KERtot(Iditi,i,ic) = kernel;
            fp(i) = mean(kernel);
            
            onexi = [ones(size(xi,1),1) xi];
            xk = onexi.*repmat(kernel, 1, p+1);
            xkonexi = xk'*onexi;
            XWY = xk' * yi;
            abi = pinv(xkonexi)*XWY;
            
            cvi = pinv(xkonexi - xk(ii,:)'*onexi(ii,:)) * ...
                (XWY - xk(ii,:)'*y(i,:));
            
            cviter(i) = real(cvi(1));
            ab(:,i) = real(abi(2:p+1));
            a(i) =  real(abi(1));
        end;
        ABtot(:,:,ic) = ab;
        Atot(:,ic) = a;
        FPtot(:,ic) = fp;
        
        CVtot(:,ic) = (cviter - y).^2;        
        qCV = quantile(CVtot(:,ic),.99);
        qChi2 = chi2inv(.99,mean(CVtot(:,ic)));
        irbCV = find(CVtot(:,ic) <= max(qCV,qChi2));
        Licv(ic) = length(irbCV);
        
        ic = ic + 1;
    end
    L0 = min(Licv);
    CVtot0 = sort(CVtot);
    CVscore = mean(CVtot0(1:L0,:));
    [~,kh] = min(CVscore);
    fp = FPtot(:,kh);
    ab = ABtot(:,:,kh);
    a = Atot(:,kh);

    normab = sum(ab.^2) + eps;
    qab1 = quantile(normab,.25);
    qab3 = quantile(normab,.75);
    qab95 = quantile(normab,.95);
    ubnd = max(qab3 + 1.5 * (qab3 - qab1),qab95);
    wab = (normab < ubnd);
    normab2 = normab .* wab + ubnd * (1 - wab);
    ab = ab .* (ones(p,1) * sqrt(normab2 ./ normab));  
    
    ay = CVtot(:,kh);
    may = mean(ay);
    qay1 = chi2inv(.99,may);
    qay2 = quantile(ay,.95);
    qay3 = quantile(ay,.975);
    qay4 = quantile(ay,.99);
    aerr = (ay <= max(qay1,qay2));
    
    fw = wt0 .* aerr;
    Ifw0 = (ay <= max(qay1,qay3));
    Ifw = (ay <= max(qay1,qay4));
    
    Mab = ab * diag(fw) *ab' / (sum(fw) + eps);
    ssMab = sum(sum(Mab));
    if isnan(ssMab) || isinf(ssMab)
        B = Btemp;
    else
        [B0 D] = eig(Mab);
        [D, I] = sort(abs(diag(D)),'descend');
        
        cD = cumsum(D / sum(D));
        nd0 = max(2,sum((cD < .9))+1);
        
        nd2 = sum(D >= D(nd0)/2);
        Itemp = I(1 : nd2);
        B = B0(:,Itemp);
        if sqrt(n) <= n * .05
            B = B .* (abs(B) > 1 / sqrt(n));
            [B,~] = qr(B,0);
        end
    end
    sizeB = size(B,2);
    WB = ones(sizeB,1);
    Bdist = norm((eye(p) - B * B') * Btemp);
    Bdist2 = norm((eye(p) - B * B') * Btemp2);
    
    if sizeB == size(Btemp,2) && sum(Bdist) < nd * .005
        break
    elseif sizeB == size(Btemp2,2) && sum(Bdist2) < nd * .0005
        Mab = B*B' + Btemp*Btemp';
        [V,D] = eig(Mab);
        [~,ID] = sort(abs(diag(D)),'descend');
        B = V(:,ID(1:sizeB));
        break
    else
        Btemp2 = Btemp;
        Btemp = B;        
    end
    if iter == 1
        Bf = B;
        af = a;
        abf = ab;
        WBf = WB;
        fpf = fp;
        Ifwf = Ifw0;
    end
end;
Ifw = Ifw0;
if iter == ITER % if fail to converge
    B = Bf;
    a = af;
    ab = abf;
    WB = WBf;
    fp = fpf;
    Ifw = Ifwf;
end

sizeB = size(B,2);
if sizeB > nd
    CV0 = 1e3*var(y);
    onen = ones(n,1);
    hp = size(B,2);
    abB = ones(hp,n);
    aB = ones(n,1);
    cviter = ones(n,1);
    XB = xx * B;
    XBWB = XB * diag(sqrt(WB));
    sigXB = min(iqr(XB)/1.34,1);
    HCV = n^(-1/(hp+4))*(2:.5:(log(n)+1)^1.5)/log(n)*log(p/2);
    for hcv = HCV
        h = hcv * sigXB;
        H = onen * h;
        spihp = (2 * pi)^(hp/2);
        
        for i = 1 : n;            
            ii = i;
            xiB = XBWB - repmat(XBWB(ii,:),n,1);
            yi = y;
            Hi = H;
            xd = xiB ./ Hi;
            
            xd2 = xd.^2;
            kernel0 = prod(exp(-xd2 / 2) ./ Hi, 2) / spihp + eps;
            kernel = kernel0;% + .5 * KerLi; %%%%%%%%%%%%%%%
            fp(i) = mean(kernel);
            
            onexi = [ones(size(xiB,1),1) xiB];
            xk = onexi.*repmat(kernel, 1, hp+1);
            xkonexi = xk'*onexi;
            XWY = xk' * yi;
            abi = pinv(xkonexi)*XWY;
            
            cvi = pinv(xkonexi - xk(ii,:)'*onexi(ii,:)) * ...
                (XWY - xk(ii,:)'*y(i,:));
            
            cviter(i) = cvi(1);
            abB(:,i) = abi(2:hp+1);
            aB(i) =  abi(1);
        end;
        
        CVtotB = (cviter - y).^2;
        qCV = quantile(CVtotB,.99);
        qChi2 = chi2inv(.99,mean(CVtotB));
        IfwB = (CVtotB <= max(qCV,qChi2));
        irbCV = find(IfwB);
        CVBscore = mean(CVtotB(irbCV,:));
        
        if CVBscore < CV0
            CV0 = CVBscore;
            abB0 = abB;
            Ifw0 = IfwB;
        end
        ic = ic + 1;
    end
    MabB = abB0 * diag(Ifw0) * abB0';
    [VB,DB] = eig(MabB);
    [~,IB] = sort(abs(diag(DB)),'descend');
    ab = B * abB0;
    B = B * VB(:,IB(1:nd));
end
sizeB = size(B,2);
WB = WB(1:sizeB);
B = real(ss*B);
for i = 1:sizeB
    B(:,i) = B(:,i)/sqrt(B(:,i)'*B(:,i) + eps);
end
ab = real(ab'*ss');