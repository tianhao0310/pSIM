function B = rMAVE(xx, y, h, nd)
%Searching the effective dimension reduction subspace of model
%  			y = g(B^TX) + e
%Useage: directions(x, y, h, d)
%Input:
%     x --- expanaltory variables
%     y --- response
%     h --- bandwidth
%     nd --- working dimension of the space
%Output:
%     Directions B
%Reference: Y. Xia, H. Tong, W.K. Li and L.X. Zhu,
% "Adaptive estimation of the effictive dimension space",(2002)
%-------------------------------------------------------------
[n,p] = size(xx);
mm = mean(xx);
xx = xx - repmat(mm,n,1);
ss = pinv(xx'*xx/n)^0.5;
xx = xx*ss;

b = (1:p)'; b = b/sqrt((b')*b);
c = b; Jzf1 = ones(n, 1); txt = Jzf1; J1 = ones(n,1); zfdd = zeros(1,p);
zfa0 = zeros(n, p*p); zfa1 = zeros(n, p); zfb0 = zeros(n, p);
zfb1 = ones(n, 1); zfb2 = ones(n, 1); nref = 10; niter = 10;
zfh0 = mean(std(xx))/n^(1/(p+4));
for zfia = 1:n;
    zfxy = xx - repmat(xx(zfia,:), n,1);
    zfx = zfxy.*zfxy*ones(p,1);
    zfK = exp(- zfx /(2*zfh0*zfh0*p) );
    K = zfK/sum(zfK); Kv = repmat(K,1,p).*zfxy;
    zfa0(zfia,:) = reshape((Kv')*zfxy, 1, p*p);
    zfa1(zfia,:) = ones(1,n,1)*Kv;
    zfb0(zfia,:) = (y')*Kv;
    zfb1(zfia) = (K')*y;
end;

zfbb0 = []; zfcc0 = zeros(p, nd);
for ik = 1 : min(nd+1,p);
    if ik > 1;
        b = (eye(p) - zfbb0*(zfbb0'))*ones(p,1);
    end;
    for iter = 1 : niter;
        b0 = b;
        zfbb = [zfbb0 b]; zfAs = 0.0; zfBs = 0.0;
        for zfia = 1 : n;
            zfaa = reshape(zfa0(zfia,:), p, p);
            zfat12 = zfa1(zfia,:)*zfbb ; zfat22 = (zfbb')*zfaa*zfbb;
            zfat1 = [1;(zfat12')];
            zfat2 = [zfa1(zfia,:)*zfbb;(zfbb')*zfaa*zfbb ];
            zfat = [zfat1 zfat2]  + (1.0E-5)*eye(1+ik);
            Szfat = sum(sum(zfat));
            if isnan(Szfat) || isinf(Szfat)
                zfbt = [mean(y);ones(ik,1) * eps];
            else
                zfbt = [zfb1(zfia)   zfb0(zfia,:) * zfbb ]';
                zfbt = pinv(zfat)*zfbt;
            end
            a = zfbt(1);
            
            zfbe1 = 0;
            if ik >1 ;
                zfbe1 = zfbt(2:ik);
            end;
            d = zfbt(1+ik);
            
            if ik == 1
                zfbb0 = zeros(p,1);
            end
            zfBs = zfBs + d*((zfb0(zfia,:)- zfa1(zfia,:)*a)' - zfaa*zfbb0*zfbe1);
            zfAs = zfAs + d*d*zfaa;
            if ik == 1
                zfbb0=[];
            end
        end;
        
        zfAs = zfAs + 1.0e-10*eye(p);
        b = pinv(zfAs)*zfBs;
        if ik > 1
            b = (eye(p) - zfbb0*(zfbb0'))*b;
        end;
        b = b/sqrt((b')*b);
        
        if  ( (b0')*(eye(p) - b*(b'))*b0 < 1.0E-5 )
            iter = niter + 1;
        end;
    end;
    if ik == 1
        zfbb0 = b;
    end
    if ik > 1
        zfbb0 = [zfbb0 b];
    end
end;
error1 = 0.0;
error2 = 0.0;
for kk = 1 : nref;
    zfbbi = zfbb0(:,1:nd); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xxk = xx*zfbb0;
    [n1,nd1] = size(xxk);
    zfh = h;
    for zfia = 1 : n;
        zfxyk = xxk - repmat(xxk(zfia,:), n,1);
        zfxy = xx - repmat(xx(zfia,:), n, 1);
        zfx = zfxyk.*zfxyk*ones(nd1,1);
        zfK = exp(- zfx /(2*zfh*zfh) );
        K = zfK/sum(zfK); Kv = repmat(K,1,p).*zfxy;
        zfa0(zfia,:) = reshape((Kv')*zfxy, 1, p*p);
        zfa1(zfia,:) = ones(1,n)*Kv;
        zfb0(zfia,:) = (y')*Kv;
        zfb1(zfia) = (K')*y;
    end;
    zfbb0 = zfbb0(:,1:nd);
    for ik = 1 : nd;
        b = zfbb0(:,ik);
        for iter = 1 : niter;
            b0 = b; zfAs = 0.0; zfBs = 0.0;
            for zfia = 1 : n;
                zfaa = reshape(zfa0(zfia,:), p, p);
                zfat12 = zfa1(zfia,:)*zfbb0; zfat22 = (zfbb0')*zfaa*zfbb0;
                
                zfat1 = [1 zfat12]';
                zfat2 = [zfa1(zfia,:)*zfbb0;(zfbb0')*zfaa*zfbb0] ;
                zfat = [zfat1 zfat2]  + (1.0E-5)*eye(1+nd);
                zfbt = [zfb1(zfia)   zfb0(zfia,:)*zfbb0]';
                
                zfbt = pinv(zfat)*zfbt;
                a = zfbt(1);
                zfbe1 = zfbt(2:(nd+1));
                d = zfbe1(ik);
                DD(zfia, ik) = d;
                zfbe1(ik) = 0.0;
                zfBs = zfBs + d*( (zfb0(zfia,:)- zfa1(zfia,:)*a)' - zfaa*zfbb0*zfbe1  );
                zfAs = zfAs + d*d*zfaa;
            end;
            zfAs = zfAs + 1.0e-10*eye(p);
            b = pinv(zfAs)*zfBs;
            zfbb0(:,ik) = zeros(p,1);
            b = (eye(p) - zfbb0*(zfbb0'))*b; b = b/sqrt((b')*b);
            zfbb0(:,ik) = b;
        end;
    end;
    error1 =  sum( diag( zfbb0'*(eye(p)-zfbbi*zfbbi')*zfbb0 ) );
    if error2 + error1 < 0.0001
        break;
    end;
    error2 = error1;
end;
B = zfbb0;

DD = DD-repmat(mean(DD), n, 1);
[v, d] = eig(DD'*DD); d = diag(d);
[d, I] = sort(d); v = v(:,I);
B = B * v;

B = real(ss * B);
B = B ./ repmat(sqrt(sum(B.^2)),p,1);
% B1 = B2;
% for i = 1:size(B,2);    
%     B1(:,i) = B1(:,i) / sqrt(B1(:,i)'*B1(:,i));
% end
% B1 = real(B1);
% B2 = real(B2);