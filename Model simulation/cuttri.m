function Tout=cuttri(T,r,L,flag)
if flag>0
    d=floor(1/(r/L));
    Tout=T;
    for n=1:r
        Tout(n,1+n*d:end)=NaN;
        Tout(end-n+1,1+n*d:end)=NaN;
    end
    Tout=Tout(:);
    Tout(isnan(Tout))=[];
elseif flag<0
    d=floor(1/(r/L));
    Tout=T;
    for n=1:r
        Tout(n,1:end-n*d)=NaN;
        Tout(end-n+1,1:end-n*d)=NaN;
    end
    Tout=Tout(:);
    Tout(isnan(Tout))=[];  
end
