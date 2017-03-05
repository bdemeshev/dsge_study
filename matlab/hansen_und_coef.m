%Method of undetermined coefficients for Hansen model
theta=0.36;
beta=0.99;
delta=0.025;
A=1.72; 
rbar=1/beta-1+delta; 
hbar=(1-theta)*(1-beta*(1-delta))/((1-theta)*(1-beta*(1-delta))+A*(1-beta*(1-delta)-beta*delta*theta));
kbar=hbar*(theta/(1/beta-1+delta))^(1/(1-theta));
ybar=kbar^theta*hbar^(1-theta);
cbar=ybar-delta*kbar; 
A=[0 -kbar 0 0]'; 
B=[0 (1-delta)*kbar theta -1]';
C=[1 -1 -1/(1-hbar) 0; ybar -cbar 0 0;-1 0 1-theta 0; 1 0 0 -1]; 
D=[0 0 1 0]' ; 
F=[0];
G=[0]; 
H=[0];
J=[0 -1 0 beta*rbar];
K=[0 1 0 0];
L=[0];
M=[0];
N=[0.95]; 
a=F-J*inv(C)*A; 
b=-(J*inv(C)*B-G+K*inv(C)*A); 
c= -K*inv(C)*B+H;
discr=b^2-4*a*c; 
p1=(-b+sqrt(discr))/(2*a);
p2=(-b-sqrt(discr))/(2*a);
P=p2
R=-inv(C)*(A*P+B)
q1=kron(N', (F-J*inv(C)*A))+(J*R+F*P+G-K*inv(C)*A); 
q2=(J*inv(C)*D-L)*N+K*inv(C)*D-M; 
Q=q2/q1 
S=-inv(C)*(A*Q+D)



