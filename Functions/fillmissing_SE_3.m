%Tommy Mitchel 2018, tmitchel@jhu.edu

function g_t = fillmissing_SE_3(g_0)

count = 0;
for i = 1:size(g_0, 3)
    
    if ~isnan(g_0(:, :, i))
        count = count + 1;
        g_exist(count) = i;
    end
end


g_t = g_0;

if count ~= 0
for i = 1:length(g_exist)-1
    if g_exist(i+1) - g_exist(i) > 1
    B = log_SE(inv_SE(g_0(:, :, g_exist(i)))*g_0(:, :, g_exist(i+1)));
    for j = g_exist(i):g_exist(i+1)
        t = ((j-g_exist(i))/(g_exist(i+1) - g_exist(i)));
        g_t(:, :, j) = g_0(:, :, g_exist(i))*exp_se(t*B);
    end
    end
end
end
        
    