% Cars
addpath('/home/svu/g0900762/pSIMRev')
clear
tic

%% Initialization
load cars
Xvar = cars(:,2 : end);
[n,p] = size(Xvar);
% X1 = Xvar(:,1);
% for i = 1:n
%     if X1(i) == 3 || X1(i) == 5
%         X1(i) = 4;
%     end
% %     X1(i) = round(X1(i)/2);
% end
% Xvar(:,1) = X1;
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

X = [Xtemp,Xdmy2];
y = cars(:,1);
Xp = [];

%% Two pieces
for m0 = 1 : 2
    m0
    [pB,yh,yp,gI,gIp,cv0,BIC,Bw] = pSIMfit(X,y,m0,Xp);
    ASE = mean((y - yh).^2)
    cv0 = cv0
    BIC
    if m0 ~= 2, continue, end
    I1 = find(gI == 1);
    I2 = find(gI == 2);
    
    XB = X * Bw;
    p1 = 1;
    p2 = 2;
    
    figure
    subplot(2,2,1)
    XpB = X * pB;
    [~,Itemp1] = sort(XpB(I1,1));
    I1 = I1(Itemp1);
    [~,Itemp2] = sort(XpB(I2,2));
    I2 = I2(Itemp2);
    plot(XpB(I1,1),y(I1,:),'o',XpB(I1,1),yh(I1,:),'r+-') % ,XpB(Ip1,1),yp(Ip1,:),'g--')
    subplot(2,2,2)
    plot(XpB(I2,2),y(I2,:),'o',XpB(I2,2),yh(I2,:),'r+-') % ,XpB(Ip2,2),yp(Ip2,:),'g--')
    subplot(2,2,3)
    plot(XB(I1,1),XB(I1,2),'o',XB(I2,1),XB(I2,2),'r+')
    subplot(2,2,4)
    Ic1 = union(find(Xvar(:,1) == 5),find(Xvar(:,1) == 6));
    Ic2 = find(Xvar(:,1) == 8);
    Ic3 = setdiff(1:n,union(Ic1,Ic2));
    plot(XB(Ic1,p1),XB(Ic1,p2),'k.',XB(Ic2,p1),XB(Ic2,p2),'r+',...
        XB(Ic3,p1),XB(Ic3,p2),'bo')
end

toc





