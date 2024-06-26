function[w,thetax,thetay]=StaticSolver(K_c,F_c,mesh,nctot)
%nctot 是约束自由度索引号
s=K_c\F_c; %删除约束自由度对应行列之后剩下的自由度的解
index=1:mesh.nn*3;
q=zeros(mesh.nn*3,1);
index(nctot)=[];%删除被约束的自由度

for i=1:numel(index)
    q(index(i))=s(i);
end
w=q(1:3:end);
thetax=q(2:3:end);
thetay=q(3:3:end);
return