function [q,S1]=struct_factor_sq(xy,Lx,Ly,L,sf,n,m)
%parameters 
% the code modifed from Mingji Huang et al PNAS 2021, 118 (18) e2100493118
% https://doi.org/10.1073/pnas.2100493118
% L=400;  % 系统大小L
% n=100;  % 长度分n份
% m=64;   % 角度分m份
% sf=0.1; % 缩放因子[0~1]



xy(:,1)=mod(xy(:,1),Lx)-Lx/2;
xy(:,2)=mod(xy(:,2),Ly)-Ly/2;
% plot(xy(:,1),xy(:,2),'.')
Lmax=L;
r1=xy;
theta=[1:m]/m*pi;       
S=[];
sec=10.^[1:log10(L)/n:log10(L)];
for ij=1:10
    if sf==0
        q=2*pi/(Lmax/sqrt(2)).*sec;
        s=zeros(length(q),length(theta));
        clear r2
        for jj=1:length(theta)
            r2(:,1)=cos(theta(jj))*r1(:,1)-sin(theta(jj))*r1(:,2); %旋转theta角后的XY坐标
            r2(:,2)=sin(theta(jj))*r1(:,1)+cos(theta(jj))*r1(:,2);
            for k=1:length(q)
                r3=r2(abs(r2(:,1))<L/2/sqrt(2)&abs(r2(:,2))<L/2/sqrt(2),:);
                s(k,jj)=abs(sum(exp(-1i*q(k)*r3(:,1)))).^2/length(r3); 
            end
        end                  
    else
        q=2*pi/(Lmax/sqrt(2)).*sec;
        s=zeros(length(q),length(theta));
        clear r2
        for jj=1:length(theta)
            r2(:,1)=cos(theta(jj))*r1(:,1)-sin(theta(jj))*r1(:,2); %旋转theta角后的XY坐标
            r2(:,2)=sin(theta(jj))*r1(:,1)+cos(theta(jj))*r1(:,2);
            for k=1:length(q)
                r3=r2(abs(r2(:,1))<L/2/sqrt(2)/(1+(ij-1)/(10./sf)) & abs(r2(:,2))<L/2/sqrt(2)/(1+(ij-1)/(10./sf)),:);
                s(k,jj)=abs(sum(exp(-1i*q(k)*r3(:,1)))).^2/length(r3); 
            end
        end
    end
    s1=nanmean(s,2);
    S=[S;s1'];
end
    S1=nanmean(S);
% figure
% loglog(q,S1,'r.')
% ylim([0.01 10])
S1=smooth(S1,5);
end

