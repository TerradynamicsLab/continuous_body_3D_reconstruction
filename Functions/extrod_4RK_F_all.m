% This function gives the rigid body motion by means of variational calculus
% of rotation and motion group in case of extensible rods.
% - made by Jin Seob Kim
% Refer to _Conformational analysis of stiff chiral polymers with
% end-constraints_ (Kim 2006) for detailed explanations
% Modified by Qiyuan to understand the codes

% Input:
%   n: should be 6 by 1 vector, xi (initial guess?)
%   L: length of segment
%   N: number of elements on the segment, recommended: >=700
%   K: stiffness matrix
%   Kinv: inverse or pseudo-inverse of K
%   kv: The k in Eqn. (11)
% Output:
%   g,ev: SE(3) matrix and twist coordinates of each element on the segment
function [g,ev]=extrod_4RK_F_all(n, L, N, K, Kinv, kv) 


%%% basis elements for SE(3)
Et1=[0,0,0,0;0,0,-1,0;0,1,0,0;0,0,0,0];
Et2=[0,0,1,0;0,0,0,0;-1,0,0,0;0,0,0,0];
Et3=[0,-1,0,0;1,0,0,0;0,0,0,0;0,0,0,0];
Et4=[0,0,0,1;0,0,0,0;0,0,0,0;0,0,0,0];
Et5=[0,0,0,0;0,0,0,1;0,0,0,0;0,0,0,0];
Et6=[0,0,0,0;0,0,0,0;0,0,0,1;0,0,0,0];

H=eye(4);
ds=L/N;
ev=zeros(6,N+1); % xi vector
ev(:,1)=n; % initial xi vector: \xi(0)=\eta
g=zeros(4*(N+1),4);
g(1:4,:)=H;

% Numerically calculate the segment using elastic rod theory
for i=2:N+1
    %%% 4th order classical Runge-Kutta method to integrate \xi
    k1=ds*fn_psi_F(ev(:,i-1), K, Kinv, kv); % fn_psi_F: Kim 2006 equation (17), d\xi/ds
    k2=ds*fn_psi_F(ev(:,i-1)+k1/2, K, Kinv, kv);
    k3=ds*fn_psi_F(ev(:,i-1)+k2/2, K, Kinv, kv);
    k4=ds*fn_psi_F(ev(:,i-1)+k3, K, Kinv, kv);
    
    ev(:,i)=ev(:,i-1)+(k1+2*k2+2*k3+k4)/6;
    
    Xmat=ev(1,i)*Et1+ev(2,i)*Et2+ev(3,i)*Et3+ev(4,i)*Et4+ev(5,i)*Et5+ev(6,i)*Et6; % Eqn. (8), g^{-1}\dot{g}
%    H=H*expm(ds*Xmat);
    H = H*exp_se(ds*Xmat);
    g(4*(i-1)+1:4*i,:)=H;
end





% %%% Penalty
% penalMat=diag([Cpen1,Cpen2,Cpen3]);
% penalvec=[0;0;0;0;0;2*Cpen3];
% 
% %%% The extensible Marko-Siggia DNA Model with penalty
% B=[a0+b0^2/c0,0,b0;0,a0,0;b0,0,c0];
% C=[0,0,0;0,0,0;0,0,tau];
% D=[s1,0,0;0,s2,0;0,0,d0];
% K=[B,C;C',D+penalMat];
% kv=K*[2*pi*nrev*[0;0;1];0;0;1]+penalvec;


% a=zeros(3,N+1);
% a(:,1)=[0;0;0];


%     h1=ds*fn_H(H,ev(:,i));
%     h2=ds*fn_H(H+h1/2,ev(:,i));
%     h3=ds*fn_H(H+h2/2,ev(:,i));
%     h4=ds*fn_H(H+h3,ev(:,i));
%     H=H+(h1+2*h2+2*h3+h4)/6;

    
%     v1=H(1:3,1); v2=H(1:3,2); 
%     v3=H(1:3,3)/norm(H(1:3,3),2);
%     v2=(v2-(v2'*v3)*v3)/norm(v2-(v2'*v3)*v3,2);
%     v1=(v1-(v1'*v3)*v3-(v1'*v2)*v2)/norm(v1-(v1'*v3)*v3-(v1'*v2)*v2,2);
%     H(1:3,1:3)=[v1,v2,v3];
%     a(:,i)=H(1:3,4);
    
