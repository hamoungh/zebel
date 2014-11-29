K=2; 
d=[1 
    1];
z=40; 
N=10; 
r=zeros(K,N); 
x=zeros(1,N);  
q=zeros(K,N+1); 
q(:,1)=[0 0]';
nn=0:9;
for n=1:10    
    r(:,n)=d.*(1+q(:,n)); 
    x(n)=n/(z+sum(r(:,n),1)); 
    q(:,n+1)=x(n) * r(:,n);
end

x
r
q