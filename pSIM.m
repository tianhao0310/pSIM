% piecewise single index model
function [pB,gI,Bpnt] = pSIM(X,y,m0,Bw,Wb,Btot,IBtot)
% Fitting the piecewise single index model
% Tested on Matlab R2010a
% Input:
%     X --- n x p matrix for expanaltory variables
%     y --- n x 1 vector for response variable
%     m0 --- number of effective pieces
%     Bw --- effective dimension reduction directions
%     Wb --- weight for Bw
%     Btot --- the local gradient directions for X
%     hOPG --- the bandwidth used to estimate Bw, Wb and Btot
% Output:
%     pB --- p x m0 matrix with each column be the estimated piecewise single-index
%     gI --- n x 1 vector for estimated group index
%     Bpnt --- n x p matrix for stable pointwise local gradient
% Reference: Wang, T. and Xia, Y,
%  "A piecewise single index model for dimension reduction",(2012)
% ---------------------------------------------------------------------- %
%% Initializing
[n,p] = size(X);
if nargin < 4
    [Bw, Wb, ~, Btot] = rOPGadap(X,y,[],min(m0,p));
end
pw = size(Bw,2);
if isempty(Wb)
    Wb = ones(pw,1);  
elseif max(Wb) < 1.5*min(Wb)
    Wb = ones(pw,1);
end

if nargin < 7
    IBtot = ones(n,1);
end
IBtot = reshape(IBtot,n,1);

mm = mean(X);
X = X - repmat(mm,n,1);
ss = real(pinv(X' * X / n)^0.5);
X = X * ss;
Bw = real(pinv(ss) * Bw);
Bw = Bw ./ repmat(sqrt(sum(Bw.^2) + eps),p,1);

Btot = real(Btot * pinv(ss));
normBtot = sqrt(sum(Btot.^2,2)) + eps;
Btot = Btot ./ repmat(normBtot,1,p);

rX = X * Bw * diag(sqrt(Wb));
knn = max([m0 * 2, fix(sqrt(n) / 5)]); 

M_Sdist = zeros(n,n);
Sdist = pdist(rX);
count = 0;
for i = 1:(n - 1)
    M_Sdist((i+1):end,i) = Sdist((1+count):(n-i+count));
    count = count + n - i;
end
M_Sdist = M_Sdist + M_Sdist';

%% Estimate local gradient direction
WtD = ones(n,1);
B = Btot;
p2 = p^2;
Bopg = zeros(n,p2);
t = fix(max([sqrt(p),log(n)]));
sqrn = sqrt(n);
logp = log(p);
flagn = (sqrt(n) <= p);
for itrn = 1 : 1
    WtD0 = WtD;
    for j = 1 : n
%         Xdist = M_Sdist(:,j);
%         [~,Ij] = sort(Xdist);
%         Ij_t = Ij(1 : t,:);
%         qntX = quantile(Xdist(Ij_t,:),.5) + eps;
%         Wtj0 = exp(-(Xdist(Ij_t,:).^2) / (2 * qntX^2));
%         Wtj = Wtj0;
%         if itrn > 1
%             SdistBj = sum((Bopg(Ij_t,:) - repmat(Bopg(j,:),t,1)).^2,2);
%             qntB = 2 * quantile(SdistBj,.5) + eps;
%             Wtj = exp(-SdistBj / qntB) .* Wtj0;
%         end
%         
%         fw = Wtj.* WtD0(Ij_t,:).*IBtot(Ij_t) + eps;
%         fw(1) = 1;
%         
%         Bj = B(Ij_t,:) .* repmat(normBtot(Ij_t,:),1,p);
%         Mtotj = Bj' * diag(fw) * Bj / (sum(fw) + eps);
%         [BMj DMj] = eig(Mtotj);
%         [DMj, Itemp1] = sort(abs(diag(DMj)));
%         B(j,:) = BMj(:,Itemp1(end))';
%         WtD(j) = DMj(end) / (sum(DMj) + eps);

%         [~,Itemp] = max(abs(B(j,:)));
        [~,Itemp] = max(abs(B(j,:)) >= 1 / sqrt(p));
        B(j,:) = B(j,:) * sign(B(j,Itemp));
        if flagn
            B(j,:) = B(j,:) .* (abs(B(j,:)) > 1 / sqrn);
        end
        B(j,:) = B(j,:) / (norm(B(j,:)) + eps);       
        
        BtB = B(j,:)' * B(j,:);
        Bopg(j,:) = reshape(BtB,1,p2);
        
        [~,Itemp] = max(abs(B(j,:)) > 1 / p);
        Btot(j,:) = Btot(j,:) * sign(Btot(j,Itemp));
    end
end
Bpnt = B;
% Bopg = B;
% [Btot,sum((Btot-Bpnt).^2,2)]
%% Cluster the points into m0 groups
IRB = find(IBtot);
BopgRB = Bopg(IRB,:);
opts = statset('MaxIter',400);
qtl = 1 / (2 * m0);
p2 = size(Bopg,2);
Cqtl = zeros(m0,p2);
for m_p = 1 : m0
%     Cqtl(m_p,:) = reshape(Bw(:,1) * Bw(:,1)',1,p2);
    Cqtl(m_p,:) = quantile(BopgRB,qtl);
    qtl = qtl + 1 / m0;
end
[gIRB,CtlRB] = kmeans(BopgRB,m0,'Options',opts,'emptyaction','drop',...
    'distance','sqEuclidean','start',Cqtl);
ng = zeros(m0,1);
for m_p = 1 : m0
    I_temp = find(gIRB == m_p);
    ng(m_p,:) = length(I_temp);
end
if min(ng) < n * min(0.15,1 / (2 * m0))
    for iter = m0 * (1 : knn) + 1
        [gIRB, CtlRB, ~, DgI] = kmeans(BopgRB,iter,'Options',opts,...
            'emptyaction','drop','distance','sqEuclidean');
        ngt = zeros(iter,1);
        for m_p = 1 : iter
            I_temp = find(gIRB == m_p);
            ngt(m_p,:) = length(I_temp);
        end
        [Sngt,Ingt] = sort(ngt,'descend');
        DgI = DgI(:,Ingt(1:m0));
        [~,gIRB] = min(DgI,[],2);
        CtlRB = CtlRB(Ingt(1:m0),:);
        if Sngt(Ingt(m0)) > n * min(0.15,1 / (2 * m0))
            break
        end
    end
end
gI = zeros(n,1);
gI(IRB) = gIRB;
NIRB = find(~IBtot);
if ~isempty(NIRB)
    for i = reshape(NIRB,1,length(NIRB))
        disti = sum((CtlRB - repmat(Bopg(i,:),m0,1)).^2,2);
        [temp,gIi] = min(disti);
        gI(i) = gIi;
    end
end

% %%%%%%%%%%%%%%%%
% gI = classify(Bopg,Bopg,gI);
% %%%%%%%%%%%%%%%%

ng = zeros(m0,1);
pB = zeros(p,m0);
for m_p = 1 : m0
    I_temp = find((gI == m_p) .* IBtot);
    ng(m_p,:) = length(I_temp);
    
    if ng(m_p,:) < log(n)
        continue;
    end
    htemp = ng(m_p,:)^(-.2);
    pB(:,m_p) = real(rMAVE(X(I_temp,:),y(I_temp,:),htemp,1));
%     pB(:,m_p) = real(rOPGadap(X(I_temp,:),y(I_temp,:),htemp,1));
    
    [~,Itemp] = max(abs(pB(:,m_p)));
    pB(:,m_p) = pB(:,m_p) * sign(pB(Itemp,m_p));
end
gI_r = gI;
% Check for swiss cheese structure
if min(ng) > log(n)
    flagS0 = zeros(n,1);
    for iter = 1 : 10
        if min(ng) < 1
            break
        end
        gI = -gI_r;
        flagS = zeros(n,1);
        for m_p = 1 : m0
            Itemp = find((gI_r == m_p));% .* IBtot);
            nI = length(Itemp);
            ng(m_p,:) = nI;
            
            Xm = X(Itemp,:) * Bw * diag(sqrt(Wb)) / sqrt(sum(Wb));
            IDXm = clusterdata(Xm,'maxclust',2);
            I1 = find(IDXm == 1); nI1 = length(I1);
            I2 = find(IDXm == 2); nI2 = length(I2);
            if isempty(I1) || isempty(I2)
                continue
            elseif min(nI1,nI2) < max(log(nI),nI / 3)
                if nI1 > nI2
                    flagS(Itemp(I2)) = 1;
                else
                    flagS(Itemp(I1)) = 1;
                end
            end
            
            Xy = [X(Itemp,:) * pB(:,m_p), y(Itemp,:)/std(y(Itemp,:))];
            IDXy = clusterdata(Xy,'maxclust',2);
            I1 = find(IDXy == 1); nI1 = length(I1);
            I2 = find(IDXy == 2); nI2 = length(I2);
            if isempty(I1) || isempty(I2)
                continue
            elseif min(nI1,nI2) < max(log(nI),nI / 3)
                if nI1 > nI2
                    flagS(Itemp(I2)) = 1;
                else
                    flagS(Itemp(I1)) = 1;
                end
            end
        end
        mng = min(ng(ng> 0));
        IflagS = find(~flagS);
        for j = 1 : n
            mg = abs(gI(j,:));
            [temp,Ij] = sort(M_Sdist(IflagS,j));
            Ij = Ij(2 : round(log(n) + knn * m0) + 1,:);
            Ij = IflagS(Ij);
            mXdist = median(M_Sdist(Ij,j));
            nj = sum(gI(Ij,:) == -mg);
            flagnj = (nj < knn * m0 * ng(mg,:) / n + 1);
            
            nk = zeros(m0,1);
            if flagnj || flagS(j)
                gItemp = gI(Ij);%%%%%%%%%%%%%%%%%%%%%%%%%
                for mgi = 1 : m0
                    Itemp = find(gItemp == -mgi);
                    mdiff = mean(sqrt(M_Sdist(Ij(Itemp),j)));
                    nwt = length(Itemp);
                    nk(mgi,:) = nwt / (mdiff + mXdist + eps) / (ng(mgi,:) + mng + eps);
                end
                [temp,gI_r(j,:)] = max(nk);
            end
        end
        if ~sum(gI + gI_r)% && iter > 1
            break
        end
        
        for m_p = 1 : m0
            I_temp = find((gI_r == m_p) .* IBtot);
            ng(m_p,:) = length(I_temp);     
            abspB = abs(pB(:,m_p));
            if p >= 20                
                Ix = find(abspB > max(abspB) / sqrt(p));%%%%%%%%%%%%%%%%%
            else
                Ix = 1 : p;%find(abspB > 1 / sqrn / 2);
            end
            
            if ng(m_p,:) < log(n)
                continue
            end
            Xm = X(I_temp,Ix);
            ym = y(I_temp,:);
            hm = ng(m_p,:)^(-.2);
            pBtemp = real(rMAVE(Xm,ym,hm,1));
            [~,Itemp] = max(abs(pBtemp));
            pBtemp = pBtemp * sign(pBtemp(Itemp));
            pB(:,m_p) = 0;
            pB(Ix,m_p) = pBtemp;
            pB(:,m_p) = pB(:,m_p) / norm(pB(:,m_p));
        end
    end
end
%% Convert pB back to the original X
Bpnt = Bpnt * ss';
pB = ss * pB;
for i = 1:size(pB,2);
    pB(:,i) = pB(:,i)/(sqrt(pB(:,i)'*pB(:,i)) + eps);
end

gI = gI_r;
for m_p = 1 : m0
    ng(m_p,:) = sum(gI == m_p);
end