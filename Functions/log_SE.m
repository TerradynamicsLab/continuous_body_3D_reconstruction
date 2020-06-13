function ln = log_SE(g)
    if size(g, 1) == 3
        R = g(1:2, 1:2);
        what = log_SO(R);
        w = so2axis(what);
        th = norm(w);
        r = g(1:2, 3);
        
        if th == 0
           ln = g - eye(3);
        else
            A = sin(th)/th;
            B = (1-cos(th))/th;
            V_inv = (1/(A^2 + B^2))*[A, B; -B, A];
            
            ln = [what, V_inv*r; zeros(1, 3)];
        end
    else
        R = g(1:3, 1:3);
        what = log_SO(R);
        w = so2axis(what);
        th = norm(w);
        
        if (th^2) == 0
            ln = g - eye(4);
        else
            r = g(1:3, 4);
            if abs(th) < 1e-3
                alpha = (1/12) + (th^2)/720 + (th^4)/30240;
            else
                A = sin(th)/th;
                B = (1 - cos(th))/(th^2);
                alpha = (1/(th^2))*(1-(A/(2*B)));
            end
            V_inv = eye(3) - (1/2)*what + alpha*(what^2);
            
            ln = [what, V_inv*r; zeros(1, 4)];
        end
    end
            