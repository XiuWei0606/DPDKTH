function D = getInnerDistance(x, L)
    dimension = length(L);
    switch dimension
        case 1
            u = L(1) - 2*x;
            l = u;
            D = l/6;
            % D = 2*l/pi^2;
        case 2
            u = L(1) - 2*x;
            v = L(2) - 2*x;
            l = 2*u*v/(u + v);
            D = l/8;
            % D = 2*l/pi^2;
        case 3
            u = L(1) - 2*x;
            v = L(2) - 2*x;
            w = L(3) - 2*x;
            l = 3*u*v*w/(u*v + v*w + u*w);
            D = l/10;
            % D = 2*l/pi^2;
    end
end