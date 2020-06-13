function Rzyx = EULERZYX_int(angles)

% angles = [alpha, beta, gamma]

Rzyx = rotz(angles(3))*roty(angles(2))*rotx(angles(1));

end 

%This rotation matrix acts on a frame by first rotating the frame aroud its principle z-axis
%by an angle gamma, then rotating the new frame around the principle
%y-axis of the original frame by an angle beta, and finally rotating
%the newest frame around the principle x-axis of the original frame by an
%angle alpha.
