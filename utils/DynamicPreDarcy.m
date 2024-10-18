function C = DynamicPreDarcy(p, dl, A, B, mu)
    nph = numel(p);
    C = cell(1, nph);
    for i = 1 : nph
        if A(i) == 0
            C{i} = 1 + 0.*p{i};
        else
            % C{i} = 1 - A(i)./(abs(p{i}) - B(i));
            C{i} = 1 - (mu{i}.*A(i))./(abs(p{i}./dl) - mu{i}.*B(i));
        end
        if i == 3
            C{i} = 1 + 0.*p{i}; % ingore gas phase
        end
    end
end