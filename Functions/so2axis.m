function w = so2axis(what)
if size(what, 1) == 2
    w = what(1, 2);
else
    w = [what(3, 2); what(1, 3); what(2, 1)];
end
