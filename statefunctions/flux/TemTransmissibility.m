classdef TemTransmissibility < StateFunction
    % Transmissibility for internal faces. May include an optional
    % pressure-dependent multiplier from a field in the fluid model.
    properties
        
    end
    
    methods
        function pp = TemTransmissibility(model)
            pp@StateFunction(model);
            if isfield(model.fluid, 'FactorNf')
                pp = pp.dependsOn('pressure', 'state');
                pp = pp.dependsOn('temperature', 'state');
                pp = pp.dependsOn({'ThermalFactorNf', 'ThermalFactorF'});
            end
            pp.label = 'T_f';
            assert(isfield(model.operators, 'T'));
            T = value(model.operators.T);
            assert(all(isfinite(T)))
            if any(T < 0)
                warning('Negative transmissibility in %d interfaces', sum(T < 0));
            end
            pp.outputRange = [0, inf];
        end

        % original
        function T = evaluateOnDomain(sfn, model, state)
            T = model.operators.T;
            if isfield(model.fluid, 'transMult')
                p = model.getProps(state, 'pressure');
                p = model.operators.faceAvg(p);
                T = model.fluid.transMult(p).*T;
            end
        end


        % test for DPDK+EDFM Thermal Permeability 
        % function T = evaluateOnDomain(prop, model, state)
        %     op = model.operators;
        %     if isfield(model.G.rock, 'THM') && model.G.rock.THM == 1
        %         %
        %         T = op.T;
        %         NumN = size(op.N, 1);
        %         indexNf = op.indexNf;
        %         indexNF = op.indexNF;
        %         indexNm = op.indexNm;
        %         flag = model.getProps(state, 'PhaseUpwindFlag');
        %         % f-f
        %         CoeNF = prop.getEvaluatedDependencies(state, 'ThermalFactorNf');
        %         % CoeNF(indexNF) = 1; CoeNF(indexNm) = 1; 
        %         % [CoeNFface, ~] = op.splitFaceCellValue(op, flag, CoeNF);
        %         CoeNFface = op.faceAvg(CoeNF);
        %         % T = T.*value(CoeNFface);
        %         % if isa(CoeNFface, 'ADI')
        %         %     T = T.*CoeNFface.val;
        %         % else
        %         %     T = T.*CoeNFface;
        %         % end
        % 
        %         % F
        %         CoeF = prop.getEvaluatedDependencies(state, 'ThermalFactorF');
        %         % Coe = CoeNF;
        %         % Coe(:) = 1;
        %         % Findex = strcmp(model.G.nnc.type, 'fracmat interior');
        %         % Findex = strcmp(model.G.nnc.type, 'fracmat boundary');
        %         % nncconection = model.G.nnc.cells(Findex, :);
        %         % Coe(nncconection(Findex, 2)) = CoeF;
        %         CoeFface = op.faceAvg(CoeF);
        %         T = T.*CoeFface;
        %         % if isa(CoeFface, 'ADI')
        %         %     T = T.*CoeFface.val;
        %         % else
        %         %     T = T.*CoeFface;
        %         % end
        % 
        %         % f-m
        %         L = op.nfSpacing; % natural fracture spacing
        %         fraction = op.fraction; % volume fractions of natural fracture and matrix
        %         y = @(x) proximity_function(x, L) - fraction(1);
        %         x0 = 0;
        %         xf = fzero(y, x0);
        %         df = xf/2;
        %         Afm = op.Afm;
        %         dm = getInnerDistance(xf, L);
        %         % permf = model.rock.perm(indexNf).*CoeNF(indexNf);
        %         permf = model.rock.perm(indexNf).*value(CoeNF(indexNf));
        %         % if isa(CoeNF, 'ADI')
        %         %     permf = model.rock.perm(indexNf).*CoeNF.val(indexNf);
        %         % else
        %         %     permf = model.rock.perm(indexNf).*CoeNF(indexNf);
        %         % end
        %         permm = model.rock.perm(indexNm);
        %         Tff = Afm.*permf./df;
        %         Tmm = Afm.*permm./dm;
        %         T((op.numNf + op.numNm + 1) : NumN) = 1./(1./Tff + 1./Tmm);
        %     else
        %         T = op.T;
        %         if isfield(model.fluid, 'transMult')
        %             p = model.getProps(state, 'pressure');
        %             p = model.operators.faceAvg(p);
        %             T = model.fluid.transMult(p).*T;
        %         end
        %     end
        % end

        % test for DP
        % function T = evaluateOnDomain(sfn, model, state)
        %     T = model.operators.T;
        %     if isfield(model.G, 'state0')
        %         op = model.operators;
        %         p = model.getProps(state, 'pressure');
        %         p_m = value(p(op.indexMapm));
        %         p_f = value(p(op.region));
        %         p_mi = model.G.state0.pressure(op.indexMapm);
        %         coefficient = (abs(p_f - p_mi) + abs(p_m - p_mi))./2./max(p_m.*0.01, abs(p_m - p_mi));
        %         % coefficient = (abs(p_f.*1.01 - p_mi) + abs(p_m.*1.01 - p_mi))./2./abs(p_m.*1.01 - p_mi);
        %         % coefficient = abs((p_f - p_mi).^2 - (p_m - p_mi).^2)./2./abs(p_f - p_m)./abs(p_m - p_mi);
        %         % coefficient = abs((p_f - p_mi).^2 - (p_m - p_mi).^2)./2./...
        %         %               max(p_m.*0.01, abs(p_f - p_m))./max(p_m.*0.01, abs(p_m - p_mi));
        %         % indm = find((p_m == p_mi));
        %         % indf = find((p_f == p_mi));
        %         % p_m(indm) = p_m(indm);
        %         % p_f(indf) = p_f(indf);
        %         % coefficient = abs(2.*p_mi - p_f - p_m)./2./abs(p_mi - p_m);
        %         % coefficient = (abs(p_f - p_mi) + abs(p_m - p_mi))./2./(abs(p_m - p_mi));
        %         % coefficient = 1.*ones(length(p_m), 1);
        %         % 1D
        %         % coefficient = (0.8374/1.6304).*abs((p_f - p_mi).^1.6304 - (p_m - p_mi).^1.6304)./ ...
        %         %               abs((p_f - p_m).*(p_m - p_mi).^(1.6304 - 1));
        %         % 2D
        %         % coefficient = (0.8394/2.0572).*abs((p_f - p_mi).^2.0572 - (p_m - p_mi).^2.0572)./ ...
        %         %               abs((p_f - p_m).*(p_m - p_mi).^(2.0572 - 1));
        %         % 3D
        %         % coefficient = (0.9322/2.2612).*abs((p_f - p_mi).^2.2612 - (p_m - p_mi).^2.2612)./ ...
        %         %               abs((p_f - p_m).*(p_m - p_mi).^(2.2612 - 1));
        %         % coefficient = (0.9322/2.2612).*abs(abs(p_f - p_mi).^2.2612 - abs(p_m - p_mi).^2.2612)./ ...
        %         %                max(p_m.*0.01, abs(p_f - p_m))./ ...
        %         %                max(p_m.*0.01, abs(p_m - p_mi)).^(2.2612 - 1);
        %         % coefficient = (0.9322/2.2612).*abs((p_f - p_mi).^2.2612 - (p_m - p_mi).^2.2612)./ ...
        %         %                max(p_m.*0.01, abs(p_m - p_mi)).^(2.2612 - 1);
        %         T(length(T) - op.numNfm + 1 : length(T)) = ...
        %             T(length(T) - op.numNfm + 1 : length(T)).*coefficient;
        %     end
        %     if isfield(model.fluid, 'transMult')
        %         p = model.getProps(state, 'pressure');
        %         p = model.operators.faceAvg(p);
        %         T = model.fluid.transMult(p).*T;
        %     end
        % end

    end
end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
