classdef MultiphaseConductiveHeatFlux < StateFunction
    %State function for computing conductive heat flux

    properties
    end

    methods
        %-----------------------------------------------------------------%
        function hf = MultiphaseConductiveHeatFlux(model, varargin)
            hf@StateFunction(model, varargin{:});
            hf = hf.dependsOn({'RockHeatTransmissibility', 'FluidHeatTransmissibility'});
            hf = hf.dependsOn({'temperature', 's'}, 'state');
            hf.label = 'H_c';
        end

        % Original
        %-----------------------------------------------------------------%
        function H = evaluateOnDomain(prop, model, state)
            op = model.operators;
            % [Tr, Tf] = prop.getEvaluatedDependencies(state, 'RockHeatTransmissibility' , ...
            %                                                 'FluidHeatTransmissibility');
            Tr = model.operators.Thr;
            Tf = model.operators.Thf;
            s = model.getProps(state, 's');
            if ~iscell(s)
                s = num2cell(s, 1);
            end
            flags = model.getProps(state, 'PhaseUpwindFlag');
            nph   = model.getNumberOfPhases();
            Tf_face = cell(1, nph);
            for i = 1:nph
                % use face saturation to make corrections for fluid heat
                % transmissbility
                Tf_face{i} = Tf(:, i).*op.faceAvg(s{i});
            end
            Tf_t = 0;
            for i = 1:nph
                % sum
                Tf_t = Tf_t + Tf_face{i};
            end
            % Tf_t = Tf_t./nph;
            T = model.getProps(state, 'temperature');
            Tt = Tr + Tf_t;
            if isfield(model.G, 'state0')
                op = model.operators;
                T_m = value(T(op.indexMapm));
                T_f = value(T(op.region));
                T_mi = model.G.state0.T(op.indexMapm);
                coefficient = (abs(T_f - T_mi) + abs(T_m - T_mi))./2./max(T_m.*0.01, abs(T_m - T_mi));
                % coefficient = (abs(T_f.*1.01 - T_mi) + abs(T_m.*1.01 - T_mi))./2./abs(T_m.*1.01 - T_mi);
                % coefficient = abs((T_f - T_mi).^2 - (T_m - T_mi).^2)./2./abs(T_f - T_m)./abs(T_m - T_mi);
                % coefficient = abs((T_f - T_mi).^2 - (T_m - T_mi).^2)./2./...
                %               max(T_m.*0.01, abs(Ts_f - T_m))./max(T_m.*0.01, abs(T_m - T_mi));
                % indm = find((T_m == T_mi));
                % indf = find((T_f == T_mi));
                % T_m(indm) = T_m(indm).*1.01;
                % T_f(indf) = T_f(indf).*1.01;
                % coefficient = abs(2*T_mi - T_m - T_f)./2./abs(T_mi - T_m);
                % coefficient = (abs(T_f - T_mi) + abs(T_m - T_mi))./2./abs(T_m - T_mi);
                % 1D
                % coefficient = (0.8374/1.6304).*abs((T_f - T_mi).^1.6304 - (T_m - T_mi).^1.6304)./ ...
                %               abs((T_f - T_m).*(T_m - T_mi).^(1.6304 - 1));
                % 2D
                % coefficient = (0.8394/2.0572).*abs((T_f - T_mi).^2.0572 - (T_m - T_mi).^2.0572)./ ...
                %               abs((T_f - T_m).*(T_m - T_mi).^(2.0572 - 1));
                % 3D
                % coefficient = (0.9322/2.2612).*abs((T_f - T_mi).^2.2612 - (T_m - T_mi).^2.2612)./ ...
                %               abs((T_f - T_m).*(T_m - T_mi).^(2.2612 - 1));
                % coefficient = (0.9322/2.2612).*abs(abs(T_f - T_mi).^2.2612 - abs(T_m - T_mi).^2.2612)./ ...
                %                max(T_m.*0.01, abs(T_f - T_m))./ ...
                %                max(T_m.*0.01, abs(T_m - T_mi)).^(2.2612 - 1);
                Tt(length(Tr) - op.numNfm + 1 : length(Tr)) = ...
                    Tt(length(Tr) - op.numNfm + 1 : length(Tr)).*coefficient;
            end
            H = -Tt.*op.Grad(T);
        end
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