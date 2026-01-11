function picks_new=Exclu_stone(picks,Sha_Index,L)
check=1;
Ttotal=Sha_Index;
picks_id=sub2ind([L,L],picks(1,:),picks(2,:));
while check==1
[~,A]=intersect(picks_id,unique(Ttotal));
        if isempty(A)==0
            picks_id(A)=ceil(rand(length(A),1)*L.^2); 
        else
            check=0;
        end
end
[picks_R,picks_C]=ind2sub([L,L],picks_id');
picks_new=[picks_R';picks_C'];
end