function [R,C]=checkstonemove(r,c,D,flag)
% flag=1 checkmin
% flag=2 checkmax
D2=zeros(length(D)+2);
D2(2:end-1,2:end-1)=D; 
D2(2:end-1,1)=D(:,end);
D2(1,2:end-1)=D(end,:);
D2(2:end-1,end)=D(:,1);
D2(end,2:end-1)=D(1,:);
D2(1,1)=D(end,end);
D2(1,end)=D(end,1);
D2(end,1)=D(1,end);
D2(end,end)=D(1,1);
r=r+1;c=c+1;
if flag==2
    errosioncheck=1;
    while errosioncheck==1
        temp=D2(r-1:r+1,c-1:c+1);
        temp=temp-temp(2,2);
        Tmax=max(max(temp));
        if Tmax>2
            id=find(temp==Tmax);
            if length(id)>=2
                A=randperm(length(id));
                id=id(A(1));
            end
            [a,b]=ind2sub([3,3],id);
%             D2(r,c)=D2(r,c)+1;
%             D2(r+(a-2),c+(b-2))=D2(r+(a-2),c+(b-2))-1;
            r=r+a-2;
            c=c+b-2;
        else
            errosioncheck=0;
            R=r-1;C=c-1;
            if R==0
               R=L;
            end
            if C=0
               C=L; 
            end
        end
    end
        
elseif flag==1
        depcheck=1;
        while depcheck==1
            temp=D2(r-1:r+1,c-1:c+1);
            temp=temp-temp(2,2);
            Tmin=min(min(temp));
            if Tmin<-2
                id=find(temp==Tmin);
                if length(id)>=2
                    A=randperm(length(id));
                    id=id(A(1));
                end
                [a,b]=ind2sub([3,3],id);
%                 D2(r,c)=D2(r,c)-1;
%                 D2(r+(a-2),c+(b-2))=D2(r+(a-2),c+(b-2))+1;
                r=r+a-2;
                c=c+b-2; 
            else
                depcheck=0;
            end
        end
end
        if c==length(D)
            D2(:,2)=D2(:,end);
        end
        if c==1
            D2(:,end-1)=D2(:,1);
        end
        if r==length(D)
            D2(2,:)=D2(end,:);
        end
        if r==1
            D2(end-1,:)=D2(1,:);   
        end
%         D_new=D2(2:end-1,2:end-1);
        
end