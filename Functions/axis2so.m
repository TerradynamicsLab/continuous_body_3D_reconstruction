function what = axis2so(w)

if length(w) == 1
    what = [0, -w; w, 0];
else
    what = [0, -w(3), w(2); w(3), 0, -w(1); -w(2), w(1), 0];
end
