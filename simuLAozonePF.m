% LA ozone
addpath('/home/svu/g0900762/pSIMRev')
clear
% initialization
tic
load LAozone
load Irnk330
Xvar = LAozone(:,2 : end);
n0 = size(Xvar,1);
Xtemp = Xvar;
mm = mean(Xtemp);
Xtemp = Xtemp - repmat(mm,n0,1);
Xtemp = Xtemp * diag(std(Xtemp).^(-1));
X0 = Xtemp;

[n0,p] = size(X0);
y0 = LAozone(:,1); 

matlabpool
ic = 1;
ITER = 100;
ncases = 1;
ASE1 = zeros(ITER,ncases);
ASEp1 = zeros(ITER,ncases);
BIC1 = zeros(ITER,ncases);
CV1 = zeros(ITER,ncases);
ASE2 = zeros(ITER,ncases);
ASEp2 = zeros(ITER,ncases);
BIC2 = zeros(ITER,ncases);
CV2 = zeros(ITER,ncases);
ASE3 = zeros(ITER,ncases);
ASEp3 = zeros(ITER,ncases);
BIC3 = zeros(ITER,ncases);
CV3 = zeros(ITER,ncases);

for np = 50
n = n0 - np;
parfor iter = 1 : ITER
    rnk = Irnk(iter,:);
    I1 = rnk(1 : n);
    I2 = rnk(n + (1 : np));    
    
    X1 = X0(I1,:);
    y1 = y0(I1,:);    

    X2 = X0(I2,:);
    y2 = y0(I2,:);
    
    % Single-index
    [pB,yh,yp,gI,gIp,cv0,BIC0,Bw,mf,fp] = pSIMfit(X1,y1,1,X2);
    ASE1(iter,ic) = mean((y1 - yh).^2);
    ASEp1(iter,ic) = mean((y2 - yp).^2);
    BIC1(iter,ic) = BIC0;

    % pSIM 2
    [pB,yh,yp,gI,gIp,cv0,BIC0,Bw,mf,fp] = pSIMfit(X1,y1,2,X2);
    ASE2(iter,ic) = mean((y1 - yh).^2);
    ASEp2(iter,ic) = mean((y2 - yp).^2);
    BIC2(iter,ic) = BIC0;
    
    % pSIM 3
    [pB,yh,yp,gI,gIp,cv0,BIC0,Bw,mf,fp] = pSIMfit(X1,y1,3,X2);
    ASE3(iter,ic) = mean((y1 - yh).^2);
    ASEp3(iter,ic) = mean((y2 - yp).^2);
    BIC3(iter,ic) = BIC0;
end
toc
end

meanASE2 = mean(ASE2)
meanASEp2 = mean(ASEp2)
medianASE2 = median(ASE2)
medianASEp2 = median(ASEp2)
save OzoneASE1 ASE1
save OzoneASEp1 ASEp1
save OzoneBIC1 BIC1
save OzoneCV1 CV1
save OzoneASE2 ASE2
save OzoneASEp2 ASEp2
save OzoneBIC2 BIC2
save OzoneCV2 CV2
save OzoneASE3 ASE3
save OzoneASEp3 ASEp3
save OzoneBIC3 BIC3
save OzoneCV3 CV3
exit