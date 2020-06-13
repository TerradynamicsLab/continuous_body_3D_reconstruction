function Tw = twist_n(w, L_seg)
N = length(w) - 1;
w_avg=sum(w(3,:))/(N+1);


Tw = (1/(2*pi))*w_avg*L_seg;



