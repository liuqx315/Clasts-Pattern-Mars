%IBM without ripple process 
%assuming that wind blows from one side

clear;
close all;
%% parameter
L=1024;     % size
sds=15;     %sands per lattice (or 25)
sal_d=2;    % saltation distance
p=0.2;      %probability to deposit
d_r=5;      %clast ratio
d_s=2*d_r+1;%clast diameter (11x11 lattices)
num_s=2000; %clast numbers
phi=2;      %maxlocal height difference
timesteps=600; %

%% initial setup
%sand
aa=0;
D=round(rand(L)).*aa+sds; % [sds,sds+1]
movedslabs=0;
    
%clast initial distribution
Spy= ceil(rand(1,num_s).*3*L/4+L/8);
Spx= ceil(rand(1,num_s).*L);
Sp=[Spy;Spx];
Sp1= ExcluBoundary(Sp',d_s/2,L);
Sp=Sp1';

%     plot(Sp(1,:),Sp(2,:),'.')
%     xlim([0,L]);ylim([0,L])
%     hold on
disp('initial done')  
 
%% main loops
% write video
v=VideoWriter('test.mp4','MPEG-4');
v.FrameRate=25;
open(v);
Results=[];

for t=1:timesteps  

% define shadow zone & foreshadow zone based on stone locations

       shadowlength=round(d_s/2/tan(15/180*pi))-d_r;  % 15 degree horizon            
       Stone_Area=Exclu_Area_id(Sp,L,d_r,0);
       Sha_Index=Exclu_Area_id(Sp,L,d_r,shadowlength);
       Sha_Area=setdiff(Sha_Index,Stone_Area);
       FSh_Index=Exclu_Area_id(Sp,L,d_r,-shadowlength);
       FSh_Area=setdiff(FSh_Index,Stone_Area);
       
% plot shadow or foreshadow area                         
%     [a,b]=ind2sub([L,L],unique(Sha_Index));
%     plot(a,b,'r.')
%     xlim([0,L]);ylim([0,L])
%     hold on
%     plot(Sp(1,:),Sp(2,:),'b.')
       
% randomly choose points   

       picks=ceil(rand(2,L^2)*L);
%
%     plot(picks(1,:),picks(2,:),'.')
%     xlim([0,L]);ylim([0,L])
%     hold on

% exclude stone area and rechoose
% picks=Exclu_stone(picks,Stone_Area,L); 
       picks=Exclu_stone(picks,Sha_Index,L); 
       R=picks(1,:);
       C=picks(2,:);
       R_s=Sp(1,:);
       C_s=Sp(2,:);
       
       
% sands movement
      for ii=1:length(R)
              RE=R(ii);
              CE=C(ii);
              if D(RE,CE)>0
                  D(RE,CE)=D(RE,CE)-1;
                  slab=1;
                  movedslabs=movedslabs+1;
              end
              
              % one sand moves
              while slab==1
                  hop=CE+sal_d;
                  if hop>L
                      hop=hop-L;
                  end
                  prob=rand();
                                    
                  %  if deposit
                  TE=sub2ind([L,L],RE,CE);
                  if  ismember(TE,Sha_Area) 
                      D(RE,hop)=D(RE,hop)+1;
                      slab=0;
                  elseif ismember(TE,FSh_Index)
                      CE=hop;
                  elseif prob<p
                      D(RE,hop)=D(RE,hop)+1;
                      slab=0;
                  else
                      CE=hop;
                  end
              end

% errosioncheck when sands moved
             D=checksandmove(L,RE,CE,D,2);
             
%  angle of deposit after moved
             RD=R(ii);
             CD=hop;
             D=checksandmove(L,RD,CD,D,1);
      end
      
%  errosioncheck of forward shadow
      for kk=1:length(FSh_Area)
            temp_id=FSh_Area(kk);
            [R1,C1]=ind2sub([L,L],temp_id);
            D=checksandmove(L,R1,C1,D,2);
      end
                
      
% sand movement caused by the maxlocal height difference
      for k=1:length(R_s)
          RE_s=R_s(k);
          CE_s=C_s(k);
              if RE_s>L-d_r
                  REp_s=RE_s+d_r-L;
              else
                  REp_s=RE_s+d_r;
              end
              if RE_s<=d_r
                  REm_s=RE_s-d_r+L;
              else
                  REm_s=RE_s-d_r;
              end
              if CE_s>L-d_r
                  CEp_s=CE_s+d_r-L;
              else
                  CEp_s=CE_s+d_r;
              end
              if CE_s<=d_r
                  CEm_s=CE_s-d_r+L;
              else
                  CEm_s=CE_s-d_r;
              end
              N=[D(REp_s,CE_s),D(REm_s,CE_s),D(RE_s,CEp_s),D(RE_s,CEm_s)]; 
              Nmin=min(N);
              id=find(N==Nmin);
                if length(id)>=2
                    A=randperm(length(id));
                    id=id(A(1));
                end
              if D(RE_s,CE_s)-Nmin > phi
                  if id==1
                      RE_s=REp_s;
                  elseif id==2
                      RE_s=REm_s;
                  elseif id==3
                      CE_s=CEp_s;
                  elseif id==4
                      CE_s=CEm_s;
                  end
              end
          stonecheck=1;
          while stonecheck==1
              if RE_s==L
                  REp_s=1;
              else
                  REp_s=RE_s+1;
              end
              if RE_s==1
                  REm_s=L;
              else
                  REm_s=RE_s-1;
              end
              if CE_s==L
                  CEp_s=1;
              else
                  CEp_s=CE_s+1;
              end
              if CE_s==1
                  CEm_s=L;
              else
                  CEm_s=CE_s-1;
              end
              
              N=[D(REp_s,CE_s),D(REm_s,CE_s),D(RE_s,CEp_s),D(RE_s,CEm_s)];
              Nmin=min(N);
              id=find(N==Nmin);
                if length(id)>=2
                    A=randperm(length(id));
                    id=id(A(1));
                end
              if -Nmin+D(RE_s,CE_s)>phi
                  if id==1
                      RE_s=REp_s;
                  elseif id==2
                      RE_s=REm_s;
                  elseif id==3
                      CE_s=CEp_s;
                  elseif id==4
                      CE_s=CEm_s;
                  end
              else
                  stonecheck=0;
              end
          end
          R_s(k)=RE_s;
          C_s(k)=CE_s;
      end
      
      Sp=[R_s;C_s];
      Sp1= ExcluBoundary(Sp',d_s/2,L);
      Sp=Sp1';

 % save data 
      Tmpresults=[Sp(2,:)' Sp(1,:)' ones(num_s,1)*t];
      Results=[Results;Tmpresults]; 

 % plot
          h=pcolor(D_plot);
          axis equal
          xlim([1,L]);ylim([1,L])
          colorbar
          set(h,'linestyle','none')
          set(gca,'Clim',[0,40])
          hold on       
          plot_stone(Sp_plot,d_r)
          hold off
          pause(0.1) 
          title(['t = ' num2str(t,'%03d') ' step'])

      f=getframe(gcf);
      writeVideo(v,f); 
        disp(t)   
end
      close(v);  

save ( ['Results_N' num2str(num_s) 'D' num2str(aa) 'R' num2str(d_s) '.mat'],'Results')