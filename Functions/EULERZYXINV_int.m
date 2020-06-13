%<3 Ratan

function RzyxInv = EULERZYXINV_int(R)

%NOTE:
%The range of this function is (-pi,pi] x [-pi/2, pi/2] x (-pi, pi], (i.e.
%the range of psi & phi is mod(2*pi) and the range of theta is mod(pi)
%since atan2(y,x) has a range of (-pi,pi] and asin(x) has a range of
%[-pi/2,pi/2].

%This function is designed under the assumption that for any rotation
%matrix R in SO(3), R can be expressed as Rzyx(alpha, beta, gamma) (i.e.
%EULERZYX([alpha, beta, gamma])) for some real angles alpha, beta, and gamma. 


%The correspondense is given by:
%R(1,1) = Rzyx(1,1) = c(beta)*c(gamma)
%R(1,2) = Rzyx(1,2) = -c(beta)*sin(gamma)
%R(1,3) = Rzyx(1,3) = s(beta)
%R(2,1) = Rzyx(2,1) = c(gamma)*s(alpha)*s(beta) + c(alpha)*s(gamma)
%R(2,2) = Rzyx(2,2) = c(alpha)*c(gamma) - s(beta)*s(gamma)*s(alpha)
%R(2,3) = Rzyx(2,3) = -c(beta)*s(alpha)
%R(3,1) = Rzyx(3,1) = s(gamma)*s(alpha) - c(gamma)*c(alpha)*s(beta) 
%R(3,2) = Rzyx(3,2) = c(gamma)*sin(alpha) + c(alpha)*s(beta)*s(gamma)
%R(3,3) = Rzyx(3,3) = c(beta)*c(alpha)

%Using the above correspondence, EULERZYXINV(R) has a unique solution (using the 
%atan2 function) for [alpha, beta, gamma] in their range described above as long 
%as R(1,3) does not equal 1 or -1. 

beta = asin(-R(3,1));

%if  abs(beta) ~= pi/2
    
    alpha = atan2(R(3,2)/cos(beta), R(3,3)/cos(beta));

    gamma = atan2(R(2,1)/cos(beta),R(1,1)/cos(beta));
    
%If R(3,1) = 1, then beta = pi/2 & R(1,1) = R(1,2) = R(2,3) = R(3,3) = 0.
% elseif beta == pi/2
%     % In this case, we find that R(2,1) = sin(gamma + alpha) & R(3,1) =
%     % -cos(gamma + alpha), so we find that gamma + alpha = atan2(R(2,1), -R(3,1)).
%     % This means that there are no unique solutions for gamma & alpha since
%     % they can be any pair of elements in (-pi, pi] x (-pi, pi] such that
%     % gamma + alpha = atan2(R(2,1), -R(3,1)). Thus there are an infinite number
%     % of possible solutions of the form [alpha, pi/2, gamma] so if we assume
%     % that alpha = 0, a possible solution is [0, pi/2, atan2(R(2,1), -R(3,1))]
%     
%     %r = rand;
%      
%     %alpha = atan2(R(2,1), -R(3,1))*(1-r);
%     
%     %gamma = atan2(R(2,1), -R(3,1))*r;
%     
%     alpha = 0;
%     
%     gamma = atan2(R(2,1),R(1,1));
%     
%     disp('No unique solution since gamma + alpha = atan2(R(2,1), -R(3,1)). One possible solution is')
%     
% %If R(3,1) = -1, then beta = -pi/2 & R(1,1) = R(1,2) = R(2,3) = R(3,3) = 0
% elseif beta == -pi/2
%     
%     % In this case, we find that R(2,1) = sin(gamma - alpha) & R(3,1) =
%     % cos(gamma - alpha), so we find that gamma - alpha = atan2(R(2,1), R(3,1)).
%     % This means that there are no unique solutions for gamma & alpha since
%     % they can be any pair of elements in (-pi, pi] x (-pi, pi] such that
%     % gamma - alpha = atan2(R(2,1), R(3,1)). Thus there are an infinite number
%     % of possible solutions of the form [alpha, -pi/2, gamma] so if we assume
%     % that alpha = 0, a possible solution is [0, -pi/2, atan2(R(2,1), R(3,1))]
%     
%    %r = rand;
%    
%    %alpha = atan2(R(2,1), R(3,1))*(1-r);
%    
%   %gamma = atan2(R(2,1), R(3,1))*r;
%    
%     alpha = 0;
%     
%     gamma = atan2(R(2,1), R(1,1));
%     
%     disp('No unique solution since gamma - alpha = atan2(R(2,1), R(3,1)).  One possible solution is')
%     
%     
% end 

RzyxInv = [alpha, beta, gamma];

end 