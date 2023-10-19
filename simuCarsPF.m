% Cars
addpath('/home/svu/g0900762/pSIMRev')
clear
% initialization
tic
load cars
load Irnk392
Xvar = cars(:,2 : end);
[n,p] = size(Xvar);
Xend = Xvar(:,end);
ndmy2 = 3;
Xdmy2 = zeros(n,ndmy2 - 1);
for i = 1 : ndmy2 - 1
    Xdmy2(:,i) = (Xend == i);
end
Xtemp = Xvar(:,1:(end-1));
mm = mean(Xtemp);
Xtemp = Xtemp - repmat(mm,n,1);
Xtemp = Xtemp * diag(std(Xtemp).^(-1));

X0 = [Xtemp,Xdmy2];
y0 = cars(:,1);

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
save CarsASE1 ASE1
save CarsASEp1 ASEp1
save CarsBIC1 BIC1
save CarsCV1 CV1
save CarsASE2 ASE2
save CarsASEp2 ASEp2
save CarsBIC2 BIC2
save CarsCV2 CV2
save CarsASE3 ASE3
save CarsASEp3 ASEp3
save CarsBIC3 BIC3
save CarsCV3 CV3
exit