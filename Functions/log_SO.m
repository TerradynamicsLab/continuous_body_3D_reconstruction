function ln = log_SO(R)

if size(R, 2) == 2
    th = atan2(R(2, 1), R(1, 1));
    ln = axis2so(th);
else
    th = acos((trace(R) - 1)/2);
    if th == 0
        ln = zeros(3, 3);
    else
        if abs(th) < 1e-3
            gamma = (1/2) + (th^2)/12 + 7*(th^4)/720;
        else
            gamma = (th/(2*sin(th)));
        end
        ln = gamma*(R - R.');
    end
end
