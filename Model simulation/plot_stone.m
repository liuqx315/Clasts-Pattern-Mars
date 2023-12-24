function plot_stone(Sp,r)
[~,b]=size(Sp);
for j=1:b
    x=Sp(2,j);
    y=Sp(1,j);
    theta=0:0.1:2*pi;
    Circle1=x+r*cos(theta);
    Circle2=y+r*sin(theta);
    plot(Circle1,Circle2,'k','linewidth',1);    
    hold on
    fill(Circle1,Circle2,'w')
end