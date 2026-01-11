function [xy1]=ExcluBoundary(xy,a,L)
    check=1;n=1;
    while check==1 && n<1e4
        dx=mod(bsxfun(@minus,xy(:,1)',xy(:,1))-L/2,L)-L/2;
        dy=mod(bsxfun(@minus,xy(:,2)',xy(:,2))-L/2,L)-L/2;
        r=abs(dx+1i*dy);
        rx=dx./r;
        ry=dy./r;
        tp=(r/2-1*a);
        tp(tp>0)=0;
        tp(isnan(tp))=0;
        tpx=tp.*rx;
        tpy=tp.*ry;
        tpx(isnan(tpx))=0;
        tpy(isnan(tpy))=0;
        xy1=round(xy+[sum(tpx,2),sum(tpy,2)]);
        xy1(xy1>L)=xy1(xy1>L)-L;
        xy1(xy1<1)=xy1(xy1<1)+L;
        if isequal(xy1,xy)
            check=0;
        else
            xy=xy1;
        end
        n=n+1;   
    end
end

