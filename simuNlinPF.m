% Nonlinear two pieces
addpath('/home/svu/g0900762/pSIMRev')
clear
% initialization
N = 100; % number of iterations
ncases = 24;
np = 50; % sample size for out-of-sample prediction

ASE1 = zeros(N,ncases);
ASE1p = zeros(N,ncases);
ASE2 = zeros(N,ncases);
ASE2p = zeros(N,ncases);
ASE3 = zeros(N,ncases);
ASE3p = zeros(N,ncases);
BIC1 = zeros(N,ncases);
BIC2 = zeros(N,ncases);
BIC3 = zeros(N,ncases);
CV1 = zeros(N,ncases);
CV2 = zeros(N,ncases);
CV3 = zeros(N,ncases);

BIC0 = zeros(N,ncases);
iBIC0 = zeros(N,ncases);
ASE0 = zeros(N,ncases);
ASEp0 = zeros(N,ncases);

% iteration
matlabpool close force
matlabpool open 8
nc = 1;
tic
for n = [200 400] % sample size
for sig = [0.01, 0.05, 0.1, 0.2] % noise level
for k = [5 10 20]; % dimension of the predictor
    p0 =  fix(sqrt(k)); % effective dimension
    p1 = fix(p0 / 2);
    be1 = zeros(k,1);
    be2 = zeros(k,1);
    be1(1 : p1) = 1 / p1;
    be2(p1 + 1 : p0) = 1 / (p0 - p1);
    B1 = (be1 + be2) / norm(be1 + be2);
    B2 = (be1 - be2) / norm(be1 - be2);
parfor count_w = 1 : N
    randn('seed',n + sig + k / 1000 + count_w)
    rand('seed',n + sig + k / 1000 + count_w)
    yyp = zeros(n + np,1);
    yyp0 = zeros(n + np,1);
    gItot = zeros(n + np,1);
    XXp = rand(n + np,k) * 2 - 1;
    for i = 1 : n + np
        Xi = XXp(i,:);
        Xi1 = Xi * be1;
        Xi2 = Xi * be2;
        if Xi2 < 0 % Xi1
            yyp0(i,:) = (Xi1 - Xi2 + 1) * exp(-(Xi1 - Xi2 + 1)^2);
            yyp(i,:) = yyp0(i,:) + sig * randn(1,1);
            gItot(i,:) = 1;
        else
            yyp0(i,:) = (Xi1 + Xi2 + 1) * exp(-(Xi1 + Xi2 + 1)^2);
            yyp(i,:) = yyp0(i,:) + sig * randn(1,1);
            gItot(i,:) = 2;
        end
    end
    X = XXp(1 : n, :);
    y = yyp(1 : n, :);
    y0 = yyp0(1 : n, :);
    gI0 = gItot(1 : n, :);
    
    Xp = XXp(n + 1 : end, :);
    yp0 = yyp0(n + 1 : end, :);
    gIp0 = gItot(n + 1 : end, :);
    
%     sig, k
%     corr(y,y0)
    
%     % Single index model fitted by rMAVE
%     hM1 = n^(-.2);
%     [yh1M,yp1M,B1M,h1M,cv1M] = rMAVEfit(X,y,1,Xp,hM1);
%     ASE1M0(count_w,nc) = mean((yh1M - y0).^2);
%     ASE1Mp(count_w,nc) = mean((yp1M - yp0).^2);
%     stdXB1M = std(X * B1M)';
%     h1M = max([h1M,n^(-.2)*stdXB1M/2],[],2);
%     BIC1(count_w,nc) = log(mean((yh1M - y).^2)) + log(n) / n / h1M * stdXB1M;
%     CV1(count_w,nc) = min(cv1M);
    
    % pSIM 1 piece (single-index model)
    m0 = 1;
    [pB,yh,yp,gI,gIp,cv0,BICm0] = pSIMfit(X,y,m0,Xp);    
    ASE1(count_w,nc) = mean((y0 - yh).^2);
    ASE1p(count_w,nc) = mean((yp0 - yp).^2);
    CV1(count_w,nc) = min(cv0);
    BIC1(count_w,nc) = BICm0;
    
    % pSIM 2 pieces
    m0 = 2;
    [pB,yh,yp,gI,gIp,cv0,BICm0] = pSIMfit(X,y,m0,Xp);    
    ASE2(count_w,nc) = mean((y0 - yh).^2);
    ASE2p(count_w,nc) = mean((yp0 - yp).^2);
    CV2(count_w,nc) = min(cv0);
    BIC2(count_w,nc) = BICm0;
    
    % pSIM 3 pieces
    m0 = 3;
    [pB,yh,yp,gI,gIp,cv0,BICm0] = pSIMfit(X,y,m0,Xp);    
    ASE3(count_w,nc) = mean((y0 - yh).^2);
    ASE3p(count_w,nc) = mean((yp0 - yp).^2);
    CV3(count_w,nc) = min(cv0);
    BIC3(count_w,nc) = BICm0;
    
    % Choose m0 based on BIC score
    [BIC0(count_w,nc), iBIC0(count_w,nc)]= ...
        min([BIC1(count_w,nc),BIC2(count_w,nc),BIC3(count_w,nc)]);    
    ASEtemp = [ASE1(count_w,nc),ASE2(count_w,nc),ASE3(count_w,nc)];
    ASEptemp = [ASE1p(count_w,nc),ASE2p(count_w,nc),ASE3p(count_w,nc)];    
    ASE0(count_w,nc) = ASEtemp(iBIC0(count_w,nc));
    ASEp0(count_w,nc) = ASEptemp(iBIC0(count_w,nc));
end
nc = nc + 1;
end
end
end
toc
save NlinASE1 ASE1
save NlinASE1p ASE1p
save NlinASE2 ASE2
save NlinASE2p ASE2p
save NlinASE3 ASE3
save NlinASE3p ASE3p
save NlinBIC1 BIC1
save NlinBIC2 BIC2
save NlinBIC3 BIC3
save NlinCV1 CV1
save NlinCV2 CV2
save NlinCV3 CV3
save NlinASE0 ASE0
save NlinASEp0 ASEp0
save NlinBIC0 BIC0
save NliniBIC0 iBIC0

exit





