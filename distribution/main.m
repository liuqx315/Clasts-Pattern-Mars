function [results]=main

%Computation of Figure 2 and Extended data Tables 3 and 4 for
%Zheng Zhu, Quan-Xing Liu, Bernard Hallet, András A. Sipos and Gábor
%Domokos:
%
%Note that sol03405_b.jpg is needed for the analysis
%
%11/16/2022
%Andras A. Sipos
%
%Implemented in Matlab R2021a
%- uses Image Processing Toolbox
%- uses Statistics and Machine Learning Toolbox
%- 

close all
figure('units','normalized','outerposition',[0 0 1 1])
set(gcf,'Color','w');
figure(1)
sp(1)=subplot(2,2,1)

%reading data and identifying the centroids of pebble shapes
I = imread('sol03405_b.jpg');       
J = im2bw(I);
K = imcomplement(J);
[label n]=bwlabel(K,8);
rp=regionprops('table',label,'Centroid');
centers = rp.Centroid;
imshow(I)
hold on
scatter(centers(:,1),centers(:,2),'r')

%Delaunay Triangulation of the measured data
DT = delaunayTriangulation(centers(:,1),centers(:,2));
triplot(DT);
axis equal
title(strcat('(a) Measured data, N=',num2str(max(size(DT.Points)))),'Interpreter','latex','Fontsize',18)
result.measured_centers=centers;

%Simulation and Delaunay Triangulation of PPP
[xx,yy]=PoissonPointProcess(5000);      %PPP with 5000 vertices
DT2 = delaunayTriangulation(xx,yy);
[xx,yy]=PoissonPointProcess(273);       %PPP with 273 vertices
DT3 = delaunayTriangulation(xx,yy);
result.PPP_centers=[xx,yy];
hold off

sp(2)=subplot(2,2,2)
scatter(xx,yy,'k')
set(gca,'XColor', 'none','YColor','none')
hold on
triplot(DT3);
axis equal
title(strcat('(b) Poisson point process, N=',num2str(max(size(DT3.Points)))),'Interpreter','latex','Fontsize',18)
hold off

%Computing the areas and edge lengthes
[areas,edgelength]=processDT(DT);
[areas2,edgelength2]=processDT(DT2);

%Kolmogorov-Smirnov tests
[h,p1,ktest1] = kstest2(areas/mean(areas),areas2/mean(areas2),'Alpha',0.05)
[h,p2,ktest2] = kstest2(edgelength/mean(edgelength),edgelength2/mean(edgelength2),'Alpha',0.05)

%Plotting the area histogram
sp(3)=subplot(2,2,3)
h1=histogram(areas/mean(areas),30,'Normalization','pdf');
h1.BinLimits=[0,7.2];
h1.NumBins=30
hold on
h2=histogram(areas2/mean(areas2),30,'Normalization','pdf');
h2.BinLimits=[0,7.2];
h2.NumBins=30

%writing result struct
result.measured_area=areas/mean(areas);
result.PPP_area=areas2/mean(areas2);
result.measured_area_histogram=h1.Values;
result.PPP_area_histogram=h2.Values;

%fitting distribution to measured areas
pd = fitdist(areas'/mean(areas),'Gamma');
ci = paramci(pd);
fitps(1,1)=ci(1,1);
fitps(1,2)=pd.a;
fitps(1,3)=ci(2,1);
fitps(1,4)=ci(1,2);
fitps(1,5)=pd.b;
fitps(1,6)=ci(2,2);
pdf1=pdf(pd,linspace(min(areas/mean(areas)),max(areas/mean(areas)),100));
plot(linspace(min(areas/mean(areas)),max(areas/mean(areas)),100),pdf1,'b','LineWidth',3)
%kernel distribution
[f,xi] = ksdensity(areas'/mean(areas)); 
plot(xi,f,'b','LineWidth',1);

%fitting distribution to PPP areas
pd = fitdist(areas2'/mean(areas2),'Gamma');
ci = paramci(pd);
fitps(2,1)=ci(1,1);
fitps(2,2)=pd.a;
fitps(2,3)=ci(2,1);
fitps(2,4)=ci(1,2);
fitps(2,5)=pd.b;
fitps(2,6)=ci(2,2);
pdf2=pdf(pd,linspace(min(areas2/mean(areas2)),max(areas2/mean(areas2)),100));
plot(linspace(min(areas2/mean(areas2)),max(areas2/mean(areas2)),100),pdf2,'r','LineWidth',3)
%kernel distribution
[f,xi] = ksdensity(areas2'/mean(areas2)); 
plot(xi,f,'r','LineWidth',1);
xlim([0,7])

legend('Measured','PP process')
title('(c) Area distribution','Interpreter','latex','Fontsize',18)
%subtitle({[strcat('Measured: std=',num2str(std(areas/mean(areas))),', min/max=',num2str(min(areas/mean(areas))/max(areas/mean(areas))))]
%    [strcat('PPP: std=',num2str(std(areas2/mean(areas2))),', min/max=',num2str(min(areas2/mean(areas2))/max(areas2/mean(areas2))))]})
hold off

%Plotting the edge length histograms
sp(4)=subplot(2,2,4)
h1=histogram(edgelength/mean(edgelength),30,'Normalization','pdf');
h1.BinLimits=[0.0,4.8];
h1.NumBins=20
hold on
h2=histogram(edgelength2/mean(edgelength2),30,'Normalization','pdf');
h2.BinLimits=[0.0,4.8];
h2.NumBins=20
result.measured_edgelength=edgelength/mean(edgelength);
result.PPP_edgelength=edgelength2/mean(edgelength2);
result.measured_edgelength_histogram=h1.Values;
result.PPP_edgelength_histogram=h2.Values;

%fitting distribution to measured edge lengths
pd = fitdist(edgelength/mean(edgelength),'Weibull')
ci = paramci(pd);
fitps(3,1)=ci(1,1);
fitps(3,2)=pd.A;
fitps(3,3)=ci(2,1);
fitps(3,4)=ci(1,2);
fitps(3,5)=pd.B;
fitps(3,6)=ci(2,2);
pdf2=pdf(pd,linspace(min(edgelength/mean(edgelength)),max(edgelength/mean(edgelength)),100));
plot(linspace(min(edgelength/mean(edgelength)),max(edgelength/mean(edgelength)),100),pdf2,'b','LineWidth',3)
%kernel distribution
[f,xi] = ksdensity(edgelength'/mean(edgelength)); 
plot(xi,f,'b','LineWidth',1);

%fitting distribution to PPP edge lengths
pd = fitdist(edgelength2/mean(edgelength2),'Weibull');
ci = paramci(pd);
fitps(4,1)=ci(1,1);
fitps(4,2)=pd.A;
fitps(4,3)=ci(2,1);
fitps(4,4)=ci(1,2);
fitps(4,5)=pd.B;
fitps(4,6)=ci(2,2);
pdf2=pdf(pd,linspace(min(edgelength2/mean(edgelength2)),max(edgelength2/mean(edgelength2)),100));
plot(linspace(min(edgelength2/mean(edgelength2)),max(edgelength2/mean(edgelength2)),100),pdf2,'r','LineWidth',3)
%kernel distribution
[f,xi] = ksdensity(edgelength2'/mean(edgelength2)); 
plot(xi,f,'r','LineWidth',1);

%writing the data into file
xlswrite('parameters.xlsx',fitps)
legend('Measured','PP process')
title('(d) Edge length distribution','Interpreter','latex','Fontsize',18)
%subtitle({[strcat('Measured: std=',num2str(std(edgelength/mean(edgelength))),', min/max=',num2str(min(edgelength/mean(edgelength))/max(edgelength/mean(edgelength))))]
%    [strcat('PPP: std=',num2str(std(edgelength2/mean(edgelength2))),', min/max=',num2str(min(edgelength2/mean(edgelength2))/max(edgelength2/mean(edgelength2))))]})
hold off

%printing the figures seperately
for i = 1:numel(sp)
    %figure('units','normalized','outerposition',[0 0 1 1])
    newfig = figure('units','normalized','outerposition',[0 0 1 1]); 
    set(gcf,'Color','w');
    axCopy = copyobj(sp(i),newfig);
    axCopy.Position = [0.13 0.11 0.775 0.815];  % default fig pos.
    set(findall(gcf,'-property','FontSize'),'FontSize',25)
end


function [xx,yy]=PoissonPointProcess(N)
%Poisson Point Process with N vertices

%Simulation window parameters
r=1;                %radius of disk
xx0=0; yy0=0;       %centre of disk

%Simulate binomial point process
pointsNumber=N; 
theta=2*pi*(rand(pointsNumber,1));  %angular coordinates
rho=r*sqrt(rand(pointsNumber,1));   %radial coordinates

%Convert from polar to Cartesian coordinates
[xx,yy]=pol2cart(theta,rho); %x/y coordinates of Poisson points

%Shift centre of disk to (xx0,yy0)
xx=xx+xx0;
yy=yy+yy0;

function [areas,edgelength]=processDT(DT)
%Processing a Delaunay Triangulatio DT

F = freeBoundary(DT)';
tri=[];
edges=[];
tri0=DT.ConnectivityList;
edges0=DT.edges;
for i=1:max(size(tri0))
    if max(ismember(tri0(i,:),F))==0
        tri=[tri;tri0(i,:)];
    end
    if max(ismember(edges0(i,:),F))==0
        edges=[edges;edges0(i,:)];
    end
end

areas = polyarea([DT.Points(tri(:,1),1),DT.Points(tri(:,2),1),DT.Points(tri(:,3),1)]',[DT.Points(tri(:,1),2),DT.Points(tri(:,2),2),DT.Points(tri(:,3),2)]');
edgelength=((DT.Points(edges(:,1),1)-DT.Points(edges(:,2),1)).^2+(DT.Points(edges(:,1),2)-DT.Points(edges(:,2),2)).^2).^(1/2);


