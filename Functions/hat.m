function N = hat(n)

a = n(1);
b = n(2);
c = n(3);

A = [0, -c, b; c, 0, -a; -b, a, 0];

N = A;

end 
