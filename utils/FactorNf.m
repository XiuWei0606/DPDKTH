% function C = FactorNf(p, T, p0, T0, Em, vm, betam, KAf, Lc, df0)
%     E = Em/(1 - 2*vm);
%     delta = abs(2.*((p - p0) - betam.*(T - T0).*E).*Lc./(E + KAf*Lc));
%     % delta = abs(2.*(- betam.*(T - T0).*E).*Lc./(E + KAf*Lc));
%     deltaMax = 1.5*df0;
%     ind = delta > deltaMax;
%     delta(ind) = deltaMax;
%     C = (df0 + delta).^3./df0.^3;
% end
% function C = FactorNf(p, T, p0, T0, Em, vm, betam, KAf, Lc, df0, indexNf, indexNm)
%     E = Em/(1 - 2*vm);
%     delta = abs(2.*((p(indexNf) - p0) - betam.*(T(indexNm) - T0).*E).*Lc./(E + KAf*Lc));
%     % delta = abs(2.*(- betam.*(T - T0).*E).*Lc./(E + KAf*Lc));
%     deltaMax = 1.5*df0;
%     ind = delta > deltaMax;
%     delta(ind) = deltaMax;
%     C = (df0 + delta).^2./df0.^2;
% end
function C = FactorNf(p, T, p0, T0, Em, vm, betam, KAf, Lc, df0, indexNf, indexNm, NumC)
    E = Em/(1 - 2*vm);
    delta = abs(2.*((p - p0) - betam.*(T - T0).*E).*Lc./(E + KAf*Lc));
    delta(indexNf) = abs(2.*((p(indexNf) - p0) - betam.*(T(indexNm) - T0).*E).*Lc./(E + KAf*Lc));
    % delta = abs(2.*(- betam.*(T - T0).*E).*Lc./(E + KAf*Lc));
    % deltaMax = 1.0*df0;
    % ind = delta > deltaMax;
    % delta(ind) = deltaMax;
    C = (df0 + delta).^3./df0.^3;
    C(~ismember((1 : NumC)', indexNf)) = 1;
end