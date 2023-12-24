clear
% the code modifed from Mingji Huang et al PNAS 2021, 118 (18) e2100493118
% https://doi.org/10.1073/pnas.2100493118
markers = {'o','s','<','d','^','>','*','p','x','h'};
COLOR=colorcube(16); COLOR(:,4)=0.5; %colorcube ;parula
FS = 18;    % Fontsize
MS = 8;     % Marksize
LW = 2;     % Linesize

select=[7 9 10];
scale=[450,450,450];% scale pix-->um

figure('position', [100 100 500 450],'color','w');
data=[];
A=[];
for num=1:length(select)
    str=['SizeXY_' num2str(select(num)) '.mat'];
    load(str);
    SizeXY=SizeXY(:,2:3)*scale(num)./1000; % pix->mm
    xy=SizeXY;

    L0=ceil(max([SizeXY(:,1) SizeXY(:,2)])); 
    L=min(L0);
    l=logspace(0.1,log10(L/2),60);
    N=100;
    var=zeros(1,length(l),N);
    rho0=size(SizeXY,1)/prod(L0);


        c=xy(:,1:2);

    for ii=1:length(l)
        for jj=1:N
    %     c=mod(bsxfun(@plus,cr,L*rand(1,2)),L);
        tp=abs(bsxfun(@minus,c,l(ii)/2+(L-l(ii))*rand(1,2)));
        var(1,ii,jj)=(sum(prod(tp<l(ii)/2,2) )/l(ii)^2-rho0)^2;
    % var(kk,ii,jj)=(sum(sum(tp.^2,2)<(l(ii)/2)^2)/(pi*l(ii)^2/4)-rho0)^2;
        end
    end
    
    V1=mean(mean(var(1,:,:),3),1)/rho0^2;
    % V2=mean(mean(var(2,:,:),3),1)/rho0^2;
    A=[A smooth(V1,5)];
    FS=16;
 
    hpf1=plot(l,smooth(V1,5), markers{num},'markersize',MS,'color',COLOR(num,:));
    hold on
    
    disp(num)
    dat=[l' smooth(V1,5)]; 
    data=[data dat];
   ylabel('$\langle\Delta n^2\rangle/n^2$','Interpreter','latex')
xlabel('length scale, $\ell$ (mm)','Interpreter','latex')
xlim([1 1000]);
xticks([1 10 100 1000]);
ylim([1e-2 2e3])
set(gca,'Position',[0.16 0.1 0.80 0.89],'XScale','log','YScale','log','linewidth',1.0,'fontsize',FS,'xminortick','on','yminortick','on',...
    'ticklength',[0.015 0.02],'TickLabelInterpreter','latex','Layer', 'top')
yy=100.*l.^-2.4;
plot(l,yy,'r-') 
%     dlmwrite(strcat(str(1:end-4),'DFALL.csv'),data,'delimiter','\t');
end

%%
% save('SizeXY_DFALLXXX.mat','data')
% load('SizeXY_DFALL.mat')
figure('position', [100 100 500 450],'color','w');
COLOR=turbo(10); %COLOR(:,4)=1.0; %colorcube ;parula; turbo
Sid=1;
for kj=1:2:17
    plot(data(:,kj),data(:,kj+1), markers{Sid},'markersize',MS,...
        'color',COLOR(Sid,:)); %'MarkerEdgeColor','black','MarkerFaceColor',
    Sid=Sid+1;
    hold on 
end

x2=logspace(0,3,100);
y2=30*x2.^(-3.0);
x3=logspace(0,3,100);
y3=20000*x3.^(-2);
plot(x3,y3,'k--','linewidth',1);
plot(x2,y2,'k--','linewidth',1);

text(8.5,0.1,'$\alpha=3.0$','Interpreter','latex','fontsize',FS,'rotation',-55);
text(200,1,'$\alpha=2.0$','Interpreter','latex','fontsize',FS,'rotation',-48);

ylabel('$\langle\Delta n^2\rangle/n^2$','Interpreter','latex')
xlabel('length scale, $\ell$ (mm)','Interpreter','latex')
xlim([1 1000]);
xticks([1 10 100 1000]);
ylim([1e-2 2e3])
set(gca,'Position',[0.13 0.11 0.83 0.87],'XScale','log','YScale','log','linewidth',1.0,'fontsize',FS,'xminortick','on','yminortick','on',...
    'ticklength',[0.015 0.02],'TickLabelInterpreter','latex','Layer', 'top')

leg1 = legend({'sol3400','sol3403\_1','sol3405\_1','sol3454','sol3403\_2','sol3404\_3','sol3404\_4','sol3405\_2','2P1354',...
    '$\ell^{-\alpha}$'},'Location','northeast','NumColumns',2);
set(leg1,'Interpreter','latex','fontsize',FS,'box','off','FontSize',FS);

save2pdf('Fig_Density_fluctuation');
