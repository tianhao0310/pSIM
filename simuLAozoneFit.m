% Mac
addpath('/home/svu/g0900762/pSIM')
clear
% initialization
tic
load LAozone
Xvar = LAozone(:,2 : end);
n0 = size(Xvar,1);
Xtemp = Xvar;
mm = mean(Xtemp);
Xtemp = Xtemp - repmat(mm,n0,1);
Xtemp = Xtemp * diag(std(Xtemp).^(-1));
X = Xtemp;
[n,p] = size(X);
y = LAozone(:,1); 

Xp = [];

for m0 = 2
    m0
    rand('seed',m0)
    % pSIM 2
%     [Bw, Wb, a, Btot] = rOPGadap(X,y,[],2);
%     Wb = sqrt(Wb); 
%     [~,~,Bw,Wb,Btot] = rOPGfit(X,y,2,[],0);
    [pB,yh,yp,gI,gIp,cv0,BIC0,Bw,mf,fp] = pSIMfit(X,y,m0,[]); %,Bw,[1,1],Btot);
    ASE = mean((y - yh).^2)
    BIC = BIC0
    cv0 = cv0
    ng = zeros(2,1);
    for m = 1 : 2
        ng(m,:) = sum(gI == m);
        if ng(m,:) == 0
            ng(m,:) = eps;
        end
    end
    if m0 ~= 2, continue, end
    
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
    plot(X(I1,7),X(I1,6),'o',X(I2,7),X(I2,6),'r+')
end
toc





