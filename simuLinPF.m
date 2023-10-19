% Comparison of psimt0
addpath('/home/svu/g0900762/pSIMRev')
clear
N = 100; % number of iterations
beta1 = [1,1,1,1,1,0,0,0,0,0]';
beta2 = [0,0,0,0,0,1,1,1,1,1]';
pB0 = zeros(10,3);
pB0(:,1) = -beta1 - sqrt(3) * beta2; pB0(:,1) = pB0(:,1) / norm(pB0(:,1));
pB0(:,2) = -beta1 + sqrt(3) * beta2; pB0(:,2) = pB0(:,2) / norm(pB0(:,2));
pB0(:,3) = beta1 / norm(beta1);

sig = 0.5 ;
k = 10;
np = 50;

ASEtot = zeros(N,3,2);
ASEptot = zeros(N,3,2);
DistBpsim = zeros(N,3,2);
BICtot = zeros(N,3,2);
CVtot = zeros(N,3,2);
Classf = zeros(1,2,2);
In = 1;
matlabpool close force
matlabpool open 8
for n = [200,400]
parfor count_w = 1 : N
    randn('seed',n + count_w)
    X0 = randn(n + np,k);
        
    y1 = -X0 * beta1 - sqrt(3) * X0 * beta2 + 1;
    y2 = -X0 * beta1 + sqrt(3) * X0 * beta2 + 1;
    y3 = 2 * X0 * beta1 + 1;
    
    c1 = (X0 * beta2 >= 0) .* ( X0 * beta2 + sqrt(3) * X0 * beta1 >= 0);
    c2 = (X0 * beta2 < 0) .* ( -X0 * beta2 + sqrt(3) * X0 * beta1 >= 0);
    
    yy0 = y1 .* c1 + y2 .* c2 + y3 .* (1 - c1 - c2);
    gI0 = c1 + c2 * 2 + (1 - c1 - c2) * 3;
    
    X = X0(1 : n, :);
    y = yy0(1 : n, :) + randn(n,1) * sig;
    y0 = yy0(1 : n, :);
    
    Xp = X0(n + 1 : n + np, :);
    yp0 = yy0(n + 1 : n + np, :);
    
    % Piecewise single index model (pSIM)
    m0 = 2;    
    [pB,yh2,yp2,gI,gIp,cv2,BIC2] = pSIMfit(X,y,m0,Xp);
    ASE2 = mean((y0 - yh2).^2);
    ASEp2 = mean((yp0 - yp2).^2);
     
    m0 = 3;    
    [pB,yh3,yp3,gI,gIp,cv3,BIC3] = pSIMfit(X,y,m0,Xp);
    ASE3 = mean((y0 - yh3).^2);
    ASEp3 = mean((yp0 - yp3).^2);
    [DistB,IB] = m2(pB0,pB);    
    DistBpsim(count_w,:,In) = vec(DistB(IB))';
     
     m0 = 4;    
    [pB,yh4,yp4,gI,gIp,cv4,BIC4] = pSIMfit(X,y,m0,Xp);
    ASE4 = mean((y0 - yh4).^2);
    ASEp4 = mean((yp0 - yp4).^2);
    
    CVtot(count_w,:,In) = [min(cv2),min(cv3),min(cv4)];
    BICtot(count_w,:,In) = [BIC2,BIC3,BIC4];
    ASEtot(count_w,:,In) = [ASE2,ASE3,ASE4];
    ASEptot(count_w,:,In) = [ASEp2,ASEp3,ASEp4];
end
In = In + 1;
end
save Lin3PASEtot3 ASEtot
save Lin3PASEptot3 ASEptot
save Lin3PDistBpsim3 DistBpsim
save Lin3PBICtot3 BICtot
save Lin3PClassf3 Classf
save Lin3PCVtot3 CVtot

exit