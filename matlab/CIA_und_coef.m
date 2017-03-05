%Method of undetermined coefficients for CIA baseline model
theta=0.36;
beta=0.99;
delta=0.025;
A=1.72; 
B=-2.5805;
gamma=0.95;
pi=0.48;
%gbar=1.02;
gbar=1;
rbar=1/beta-1+delta; 
wbar=(1-theta)*(rbar/theta)^(theta/(theta-1)); 
cbar=-beta*wbar/gbar/B;
kbar=cbar/(rbar/theta-delta);
hbar=(rbar/theta)^(1/(1-theta))*kbar;
pbar=1/cbar;
k=2;

A=[kbar 0  0 0]'; 
B=[-(rbar+1-delta)*kbar 1-theta -theta 0]';
C=[-rbar*kbar  -wbar*hbar -wbar*hbar -1/pbar; 1 0 theta-1 0;0 1 theta 0; 0 -1 0 -1]; 
D=[0 0;-1 0;-1 0;0 pi] ; 
F=[0];
G=[0]; 
H=[0];
J=[beta*rbar -1 0 0];
K=[0 1 0 0];
L=[0 0];
M=[0 0];
N=[gamma 0;0 pi]; 
a=F-J*inv(C)*A; 
b=-(J*inv(C)*B-G+K*inv(C)*A); 
c= -K*inv(C)*B+H;
discr=b^2-4*a*c; 
p1=(-b+sqrt(discr))/(2*a)
p2=(-b-sqrt(discr))/(2*a)
if abs(p1)<1 
    P=p1;
elseif abs(p2)<1 
    P=p2;
else
    display('both roots are unstable')
    %break
end 
R=-inv(C)*(A*P+B)
q1=kron(N', (F-J*inv(C)*A))+kron(eye(k),(J*R+F*P+G-K*inv(C)*A)); 
[r,c]=size((J*inv(C)*D-L)*N+K*inv(C)*D-M);
q2=reshape(((J*inv(C)*D-L)*N+K*inv(C)*D-M),c,r); 
Q=inv(q1)*q2;
Q=Q'
S=-inv(C)*(A*Q+D)



