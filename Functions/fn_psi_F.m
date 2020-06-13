function evdot=fn_psi_F(ev, K, Kinv, kv) %,K,Kinv,kv

evdot=-Kinv*wedge(K*ev-kv,ev);






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

% if rcond(K)>1e-12

%     evdot=-inv(K)*wedge(K*ev-kv,ev);
%     evdot=-inv(K)*Liealg(K*ev-kv,ev); % same as "wedge"
%     evdot=-inv(K)*infrigidmotion(ev,K*ev-kv);
% else
%     evdot=Kinv*(ad_se3(matr(ev))'*(K*ev-kv));
%     evdot=-pinv(K)*wedge(K*ev-kv,ev);
%     evdot=-pinv(K)*Liealg(K*ev-kv,ev); % same as "wedge"
%     evdot=-pinv(K)*infrigidmotion(ev,K*ev-kv);
% end
