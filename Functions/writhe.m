%Original code by Jin Seob Kim 2006
%Modifications by Thomas Mitchel 2017

function Wr = writhe(g_seg, L_seg)

L = L_seg;
N = (length(g_seg)/4)-1;


ds = L/N;
%===== construct r and v =====%
r=zeros(3,N); % position vector
v=zeros(3,N); % tangent vector OR velocity vector

for i=1:N
    r(:,i)=g_seg(4*(i-1)+1:4*i-1,4);
    v(:,i)=g_seg(4*(i-1)+1:4*i-1,3);
end

%========== calculating writhe ==========%
sum1=0;
for i = 2:N+1
    sum2=0;
    for j = i+1:N
       sum2=sum2+ds*cross(v(:,i),v(:,j))'*(r(:,i)-r(:,j))/norm(r(:,i)-r(:,j),2)^3; 
       
    end
     sum1=sum1+ds*sum2;
end

Wr = 2*sum1/(4*pi);

