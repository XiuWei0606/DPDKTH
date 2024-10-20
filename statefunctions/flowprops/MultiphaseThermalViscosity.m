classdef MultiphaseThermalViscosity < StateFunction
%State function for computing pressure- and temperature-dependent viscosity

    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = MultiphaseThermalViscosity(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'pressure'}, 'state');
            gp = gp.dependsOn({'temperature'}, 'state');
            % if model.getNumberOfPhases() == 2
            %     gp = gp.dependsOn({'pressure', 'enthalpy'}, 'state');
            % end
            gp.label = '\mu_\alpha';
        end
        
        %-----------------------------------------------------------------%
        function mu = evaluateOnDomain(prop, model, state)
            % Get phase pressures, mass fractions and temperature
            % p = prop.getEvaluatedDependencies(state, 'PhasePressures');
            p = model.getProp(state, 'pressure');
            T = model.getProp(state, 'temperature');
            phases = model.getPhaseNames();
            mu     = cell(1, numel(phases));
            % Compute density using fluid function handle
            for i = 1:numel(phases)
                % mu{i} = prop.evaluateFluid(model, ['mu', phases(i)], p{i}, T);
                mu{i} = prop.evaluateFluid(model, ['mu', phases(i)], p, T);
            end
        end
        
        %-----------------------------------------------------------------%
        % function mu = computePhaseViscosities(prop, model, state)
        %     % Get phase pressures, mass fractions and temperature
        %     [p, X] = prop.getEvaluatedDependencies(state, 'PhasePressures'             , ...
        %                                                   'ComponentPhaseMassFractions');
        %     T = model.getProp(state, 'temperature');
        %     % Identify NaCl (if present)
        %     cnames = model.getComponentNames();
        %     ix     = strcmpi(cnames, 'NaCl');
        %     if any(ix), X = X{ix}; else, X = []; end
        %     phases = model.getPhaseNames();
        %     mu     = cell(1, numel(phases));
        %     % Compute density using fluid function handle
        %     for i = 1:numel(phases)
        %         mu{i} = prop.evaluateFluid(model, ['mu', phases(i)], p{i}, T, X);
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