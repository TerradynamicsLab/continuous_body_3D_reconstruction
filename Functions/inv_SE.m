function g_verse = inv_SE(g)

g_verse = [g(1:3, 1:3)', -(g(1:3, 1:3)')*g(1:3, 4); zeros(1, 3), 1];

