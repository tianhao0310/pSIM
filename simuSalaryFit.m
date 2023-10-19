% Hitters salary data
addpath('/home/svu/g0900762/pSIMRev')
clear
% initialization
tic
load salaryRemo3
X = salaryRemo(:,1:16);
y = log(salaryRemo(:,17));
% load salary
% X = salary(:,1:16);
% y = log(salary(:,17));

[n,p] = size(X);
Xtemp = X;
mm = mean(Xtemp);
Xtemp = Xtemp - repmat(mm,n,1);
Xtemp = Xtemp * diag(std(Xtemp).^(-1));
X = Xtemp;

Xp = [];

for m0 = 1:5 % number of pieces
    m0
    rand('seed',m0)
    [pB,yh,yp,gI,gIp,cv0,BIC0,Bw] = pSIMfit(X,y,m0,[]);
    ASE = mean((y - yh).^2)
    BIC = BIC0
    cv0 = cv0
    if m0 ~= 2, continue, end
    ng = zeros(2,1);
    IpB = [2,1];
    pB = pB(:,IpB);
    gI = 3 - gI;
    for m = 1 : 2
        ng(m,:) = sum(gI == m);
        if ng(m,:) == 0
            ng(m,:) = eps;
        end
    end
    
    pB,Bw
    I1 = find(gI == 1);
    I2 = find(gI == 2);
    XpB = X * pB;
    [~,Itemp1] = sort(XpB(I1,1));
    I1 = I1(Itemp1);
    [~,Itemp2] = sort(XpB(I2,2));
    I2 = I2(Itemp2);
    XB = -X * Bw;
    figure
    subplot(2,2,1)
    plot(XpB(I1,1),y(I1,:),'o',XpB(I1,1),yh(I1,:),'r+-')
    subplot(2,2,2)
    plot(XpB(I2,2),y(I2,:),'o',XpB(I2,2),yh(I2,:),'r+-')
    subplot(2,2,3)    
    plot(XB(I1,1),XB(I1,2),'o',XB(I2,1),XB(I2,2),'r+')
    subplot(2,2,4)
    plot(X(I1,7),y(I1,1),'o',X(I2,7),y(I2,1),'r+')
end
toc





