function g = exp_se(Xihat)

if size(Xihat, 2) == 3
    th = so2axis(Xihat(1:2, 1:2));
    if th == 0
        g = eye(3) + Xihat;
    else
        v = Xihat(1:2, 3);
        V = (1/th)*[sin(th), cos(th) - 1; 1 - cos(th), sin(th)];
        
        g = [exp_so(Xihat(1:2, 1:2)), V*v; zeros(1, 2), 1];
    end
else
    what = Xihat(1:3, 1:3);
    w = so2axis(what);
    th = norm(w);
    
    if (th^2) == 0
        g = eye(4) + Xihat;
    else
        v = Xihat(1:3, 4);
        if abs(th) < 1e-3
        A = 1 - (th^2)/6 + (th^4)/120;
        B = (1/2) - (th^2)/24 + (th^4)/720;
        C = (1/6) - (th^2)/120 + (th^4)/5040;
        else
        A = sin(th)/th;
        B = (1 - cos(th))/(th^2);
        C = (1-A)/(th^2);
        end
        R = eye(3) + A*what + B*(what^2);
        V = eye(3) + B*what + C*(what^2);
        
        g = [R, V*v; zeros(1, 3), 1];
    end
end

    
    
        
        
    