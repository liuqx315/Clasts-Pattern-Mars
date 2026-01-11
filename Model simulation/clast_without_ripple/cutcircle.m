function Tout=cutcircle(T,C,r)
[a,b]=size(T);
Tout=T;
for nx=1:a
    for ny=1:b
        if abs((nx-C)+1i*(ny-C))<r+1
            Tout(nx,ny)=T(nx,ny);
        else
             Tout(nx,ny)=NaN;
        end
    end
end
Tout=Tout(:);
Tout(isnan(Tout))=[];
end
