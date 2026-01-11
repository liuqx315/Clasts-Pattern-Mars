function Ttotal=Exclu_Area_id(Sp,L,d_r,shadowlength)

if shadowlength>0
    r=Sp(1,:);c=Sp(2,:);
    Sp_id=sub2ind([L,L],r,c);
    Ttotal=[];
    for i=1:length(Sp_id)        
        temp=ones(2*d_r+1,2*d_r+1+shadowlength).*Sp_id(i);
        c1=-d_r:1:d_r;
        d1=-d_r*L:L:(d_r+shadowlength)*L;
        T1=bsxfun(@plus,temp,c1');
            if r(i)<= d_r
                d=1:(d_r-(r(i)-1));
                T1(:,d)=T1(:,d)+L;
            end
            if r(i)>=L-d_r
                d=(d_r+L-r(i)+2):2*d_r+1;
                T1(:,d)=T1(:,d)-L;
            end
        % boundry
        T2=bsxfun(@plus,T1,d1);
        T2(T2<1)=T2(T2<1)+L^2;
        T2(T2>L^2)=T2(T2>L^2)-L^2;
        Tout=cuttri(T2(:,d_r+1:end),d_r,d_r+shadowlength,1);
        Ttotal=[Ttotal;Tout];
    end
elseif shadowlength<0
    r=Sp(1,:);c=Sp(2,:);
    Sp_id=sub2ind([L,L],r,c);
    Ttotal=[];
    for i=1:length(Sp_id)        
        temp=ones(2*d_r+1,2*d_r+1+abs(shadowlength)).*Sp_id(i);
        c1=-d_r:1:d_r;
        d1=(-d_r+shadowlength)*L:L:d_r*L;
        T1=bsxfun(@plus,temp,c1');
            if r(i)<= d_r
                d=1:(d_r-(r(i)-1));
                T1(:,d)=T1(:,d)+L;
            end
            if r(i)>=L-d_r
                d=(d_r+L-r(i)+2):2*d_r+1;
                T1(:,d)=T1(:,d)-L;
            end
        % boundry
        T2=bsxfun(@plus,T1,d1);
        T2(T2<1)=T2(T2<1)+L^2;
        T2(T2>L^2)=T2(T2>L^2)-L^2;
        Tout=cuttri(T2(:,1:end-d_r),d_r,d_r-shadowlength,-1);
        Ttotal=[Ttotal;Tout];
    end 
elseif shadowlength==0
    r=Sp(1,:);c=Sp(2,:);
    Sp_id=sub2ind([L,L],r,c);
    Ttotal=[];
    for i=1:length(Sp_id)        
        temp=ones(2*d_r+1,2*d_r+1+shadowlength).*Sp_id(i);
        c1=-d_r:1:d_r;
        d1=-d_r*L:L:(d_r+shadowlength)*L;
        T1=bsxfun(@plus,temp,c1');
            if r(i)<= d_r
                d=1:(d_r-(r(i)-1));
                T1(:,d)=T1(:,d)+L;
            end
            if r(i)>=L-d_r
                d=(d_r+L-r(i)+2):2*d_r+1;
                T1(:,d)=T1(:,d)-L;
            end
        % boundry
        T2=bsxfun(@plus,T1,d1);
        T2(T2<1)=T2(T2<1)+L^2;
        T2(T2>L^2)=T2(T2>L^2)-L^2;
        Tout=cutcircle(T2,d_r+1,d_r);
        Ttotal=[Ttotal;Tout];
    end
end

