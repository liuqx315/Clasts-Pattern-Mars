%IBM 
%assuming that wind blows from the left
%rules:  
%
function main
%% parameter
L=512;      %system size
sds=15;     %sands per lattice 
sal_d=3;    % saltation distance
p=0.2;      %probability to deposit
d_r=5;      %stone ratio
d_s=2*d_r+1;%stone diameter (11x11 lattices)
num_s=100;   %stone numbers
phi=2;      %maxlocal height difference
bedmove=1;  %consider the bed movement = 1, else= 0
timesteps=1000;
%% initial setup
%sand
aa=0;
D=round(rand(L)).*aa+sds; % [sds,sds+1]
movedslabs=0;
%DD=zeros(L,L,timesteps);    
%stone
Sp= ceil(rand(2,num_s).*L/4+3*L/8);
Sp1= ExcluBoundary(Sp',d_s/2,L);
Sp=Sp1';
%     plot(Sp(1,:),Sp(2,:),'.')
%     xlim([0,L]);ylim([0,L])
%     hold on
disp('initial done')   

%% main loops
v=VideoWriter('test.mp4','MPEG-4');
v.FrameRate=25;
open(v);
Results=[];

for t=1:timesteps 
    ka=t/timesteps;
    tic

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
%     plot(picks(1,:),picks(2,:),'.')
%     xlim([0,L]);ylim([0,L])
%     hold on

% exclude stone area and rechoose
%      picks=Exclu_stone(picks,Stone_Area,L); 
       picks=Exclu_stone(picks,Sha_Index,L); 
       R=picks(1,:);
       C=picks(2,:);
       R_s=Sp(1,:);
       C_s=Sp(2,:);
       

       
% sands move  
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
                  
               % define shadow area of bed movement
               if bedmove==1
                    shadow1=hop-1;
                    shadow2=hop-2;
                    shadow3=hop-3;
                    shadow4=hop-4;
                    shadow1(shadow1<1)=shadow1(shadow1<1)+L;
                    shadow2(shadow2<1)=shadow2(shadow2<1)+L;
                    shadow3(shadow3<1)=shadow3(shadow3<1)+L;
                    shadow4(shadow4<1)=shadow4(shadow4<1)+L;
               end
                                       
                  %  if deposit
                  TE=sub2ind([L,L],RE,CE);
                  if  ismember(TE,Sha_Area) 
                      D(RE,hop)=D(RE,hop)+1;
                      slab=0;
                  elseif ismember(TE,FSh_Index)
                      CE=hop;
                  elseif bedmove==1 && (D(RE,shadow1)-D(RE,hop)>=1 ||  D(RE,shadow2)-D(RE,hop)>=2 ...
                          ||  D(RE,shadow3)-D(RE,hop)>=3 ||  D(RE,shadow4)-D(RE,hop)>=4)
                      D(RE,hop)=D(RE,hop)+1;
                      slab=0;
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
                
      
% stone movement
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
      aa=(1+rand(num_s,1).*ka/2).*d_s/2;
      Sp1= ExcluBoundary(Sp',0.5*(aa+aa'),L);
      Sp=Sp1';
 % save data 
      Tmpresults=[Sp(2,:)' Sp(1,:)' ones(num_s,1)*t];
      Results=[Results;Tmpresults]; 
 % plot

          D_plot=D;
          Sp_plot=Sp;

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
          a=round(L/4/100)*100;
          xticks([1 a:a:L])
          yticks([1 a:a:L])
          set(gca,'Tickdir','out','Ticklength',[0.03 0.03])

        f=getframe(gcf);
        writeVideo(v,f); 
        disp(t)
        toc
%         DD(:,:,t)=D;
%         save ( ['Results_' num2str(bedmove) 'N' num2str(num_s) 'L' num2str(L) 'R' num2str(d_s) '.mat'],'Results','DD')
end
    close(v);  
end
