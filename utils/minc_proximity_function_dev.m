function dp = minc_proximity_function_dev(x, L)
    dimension = length(L);
    switch dimension
        case 1
            dp = 2/L(1);
        case 2
            dp = (-8*x + 2*L(1) + 2*L(2))/(L(1)*L(2));
        case 3
            dp = (24*x^2 - 8*(L(1)+L(2)+L(3))*x +...
                2*(L(1)*L(2) + L(1)*L(3) + L(2)*L(3)))...
                /(L(1)*L(2)*L(3));
        otherwise
            error('Too many entries for fracture spacing L');
    end
end
