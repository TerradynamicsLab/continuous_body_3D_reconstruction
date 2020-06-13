function R = exp_so(what)
if size(what, 1) == 2
    th = what(2, 1);
    R = [cos(th), -sin(th); sin(th), cos(th)];
else
    th = norm(so2axis(what));
    
    if (th^2) == 0
        R = eye(3);
    else
        if abs(th) < 1e-3
             A = 1 - (th^2)/6 + (th^4)/120;
             B = (1/2) - (th^2)/24 + (th^4)/720;
        else
             A = sin(th)/th;
             B = (1 - cos(th))/(th^2);
        end
     R = eye(3) + A*what + B*(what^2);
    end
end

