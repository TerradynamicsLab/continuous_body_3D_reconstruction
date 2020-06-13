function z=wedge(x,y)
% x, y both should be 6 x 1 vectors.
% - written by Jin Seob Kim

w1=x(1:3);
w2=y(1:3);
v1=x(4:6);
v2=y(4:6);

z=[cross(w2,w1)+cross(v2,v1);cross(w2,v1)];
