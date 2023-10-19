% Hitters' salary
addpath('/home/svu/g0900762/pSIMRev')
clear
% initialization
tic
load salaryRemo3
load Irnk258
X0 = salaryRemo(:,1:16);
y0 = log(salaryRemo(:,17));
[n0,p] = size(X0);

Xtemp = X0;
mm = mean(Xtemp);
Xtemp = Xtemp - repmat(mm,n0,1);
Xtemp = Xtemp * diag(std(Xtemp).^(-1));
X0 = Xtemp;

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

CV = zeros(ITER,ncases);
for np = 50
n = n0 - np;
parfor iter = 1 : ITER          
    rnk = Irnk(iter,:);
    Ipred = rnk(n + (1 : np));
    Ifit = rnk(1 : n);
    
    X1 = X0(Ifit,:);
    y1 = y0(Ifit,:);  
    
    X2 = X0(Ipred,:);
    y2 = y0(Ipred,:);
    
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
ic = ic + 1;
end

meanASE2 = mean(ASE2)
meanASEp2 = mean(ASEp2)
medianASE2 = median(ASE2)
medianASEp2 = median(ASEp2)
save SalaryASE1 ASE1
save SalaryASEp1 ASEp1
save SalaryBIC1 BIC1
save SalaryCV1 CV1
save SalaryASE2 ASE2
save SalaryASEp2 ASEp2
save SalaryBIC2 BIC2
save SalaryCV2 CV2
save SalaryASE3 ASE3
save SalaryASEp3 ASEp3
save SalaryBIC3 BIC3
save SalaryCV3 CV3

exit





