function D_new=checksandmove(L,r,c,D,flag)
% flag=1 checkmin
% flag=2 checkmax
if flag==2
    RE=r;CE=c;
    errosioncheck=1;
    while errosioncheck==1
        if RE==L
            REp=1;
        else
            REp=RE+1;
        end
        if RE==1
            REm=L;
        else
            REm=RE-1;
        end
        if CE==L
            CEp=1;
        else
            CEp=CE+1;
        end
        if CE==1
            CEm=L;
        else
            CEm=CE-1;
        end
        
        N=[D(REp,CE),D(REm,CE),D(RE,CEp),D(RE,CEm)];
        Nmax=max(N);
        id=find(N==Nmax);
        if length(id)>=2
            A=randperm(length(id));
            id=id(A(1));
        end
        if Nmax-D(RE,CE)>2
            D(RE,CE)=D(RE,CE)+1;
            if id==1
                D(REp,CE)=D(REp,CE)-1;
                RE=REp;
            elseif id==2
                D(REm,CE)=D(REm,CE)-1;
                RE=REm;
            elseif id==3
                D(RE,CEp)=D(RE,CEp)-1;
                CE=CEp;
            elseif id==4
                D(RE,CEm)=D(RE,CEm)-1;
                CE=CEm;
            end
        else
            errosioncheck=0;
        end
    end
    D_new=D;
    
elseif flag==1
    depcheck=1;
    RD=r;CD=c;
    while depcheck==1
        if RD==L
            RDp=1;
        else
            RDp=RD+1;
        end
        if RD==1
            RDm=L;
        else
            RDm=RD-1;
        end
        if CD==L
            CDp=1;
        else
            CDp=CD+1;
        end
        if CD==1
            CDm=L;
        else
            CDm=CD-1;
        end
        
        N=[D(RDp,CD),D(RDm,CD),D(RD,CDp),D(RD,CDm)];
        Nmin=min(N);
        id=find(N==Nmin);
        if length(id)>=2
            A=randperm(length(id));
            id=id(A(1));
        end
        if -Nmin+D(RD,CD)>2
            D(RD,CD)=D(RD,CD)-1;
            if id==1
                D(RDp,CD)=D(RDp,CD)+1;
                RD=RDp;
            elseif id==2
                D(RDm,CD)=D(RDm,CD)+1;
                RD=RDm;
            elseif id==3
                D(RD,CDp)=D(RD,CDp)+1;
                CD=CDp;
            elseif id==4
                D(RD,CDm)=D(RD,CDm)+1;
                CD=CDm;
            end
        else
            depcheck=0;
        end
    end
    D_new=D;
end
end
              