function R = axis2rot(w)
R = exp_so(axis2so(w));
