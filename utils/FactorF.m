function C = FactorF(p, T, p0, T0, Em, vm, betam, KAF, LcF, dF0, indexNF, indexNm, NumC)
    E = Em/(1 - 2*vm);
    % delta = abs(2.*((p - p0) - betam.*(T - T0).*E).*LcF./(E + KAF*LcF));
    delta = abs(2.*( - betam.*(T - T0).*E).*LcF./(E + KAF*LcF));
    % delta(indexNF) = abs(2.*((p(indexNF) - p0) - betam.*(T(indexNm) - T0).*E).*LcF./(E + KAF*LcF));
    delta(indexNF) = abs(2.*( - betam.*(T(indexNm) - T0).*E).*LcF./(E + KAF*LcF));
    % deltaMax = 1.0*dF0;
    % ind = delta > deltaMax;
    % delta(ind) = deltaMax;
    C = (dF0 + delta).^3./dF0.^3;
    C(~ismember((1 : NumC)', indexNF)) = 1;
end