function[w,thetax,thetay]=StaticSolver(K_c,F_c,mesh,nctot)
%nctot ��Լ�����ɶ�������
s=K_c\F_c; %ɾ��Լ�����ɶȶ�Ӧ����֮��ʣ�µ����ɶȵĽ�
index=1:mesh.nn*3;
q=zeros(mesh.nn*3,1);
index(nctot)=[];%ɾ����Լ�������ɶ�

for i=1:numel(index)
    q(index(i))=s(i);
end
w=q(1:3:end);
thetax=q(2:3:end);
thetay=q(3:3:end);
return