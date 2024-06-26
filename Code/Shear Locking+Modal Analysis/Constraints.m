function[K_c,M_c,F_c,nctot]=Constraints(nc,K,M,F)
nctot=zeros(1,numel(nc)*3);
for i=1:numel(nc)
    nctot(i*3-2:i*3)=[nc(i)*3-2 nc(i)*3-2 nc(i)*3-2];
end
K(nctot,:)=[];
K(:,nctot)=[];
K_c=K;
M(nctot,:)=[];
M(:,nctot)=[];
M_c=M;
F(nctot)=[];
F_c=F;
return