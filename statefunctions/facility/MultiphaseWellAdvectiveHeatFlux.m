classdef MultiphaseWellAdvectiveHeatFlux < StateFunction
%State function for advective heat flux between wellbore and reservoir
    
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = MultiphaseWellAdvectiveHeatFlux(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'PhaseFlux', 'FacilityWellMapping'});
            gp.label = 'q_{avd}';
        end
        
        %-----------------------------------------------------------------%
        function q = evaluateOnDomain(prop, model, state)
            % Get well phase fluxes and well mapping
            [v, map] = prop.getEvaluatedDependencies(state, 'PhaseFlux', 'FacilityWellMapping');
            % Get enthalpy and density in perforated cells
            % [h, rho] = model.ReservoirModel.getProps(state, 'PhaseEnthalpy', 'Density');
            % h   = cellfun(@(h)   h(map.cells)  , h  , 'UniformOutput', false);
            % rho = cellfun(@(rho) rho(map.cells), rho, 'UniformOutput', false);
            % Using PhaseInternalEnergy
            [u, rho] = model.ReservoirModel.getProps(state, 'PhaseInternalEnergy', 'Density');
            u   = cellfun(@(u)   u(map.cells)  , u  , 'UniformOutput', false);
            rho = cellfun(@(rho) rho(map.cells), rho, 'UniformOutput', false);
            % Identify injecting perforations
            vt = 0;
            for i = 1:numel(v)
                vt = vt + v{i};
            end
            injector = value(vt) > 0;
            if any(injector)
                % For injecting perforations, compute enthalpy as seen from
                % wellore and update advecitve heat flux accordingly
                % Get bottom-hole pressure from well variables
                isbhp = strcmpi(state.FacilityState.names, 'bhp');
                if any(isbhp)
                    bhp = state.FacilityState.primaryVariables{isbhp}(map.perf2well);
                else
                    bhp = vertcat(state.wellSol(map.active).bhp);
                    bhp = bhp(map.perf2well);
                end
                % Account for gravity
                cp  = bhp + vertcat(state.wellSol(map.active).cdp);
                % Get temperature from well variables
                isT = strcmpi(state.FacilityState.names, 'well_temperature');
                if any(isT)
                    Tw = state.FacilityState.primaryVariables{isT}(map.perf2well);
                else
                    Tw = vertcat(map.W.T);
                    Tw = Tw(map.perf2well);
                end
                % Compute enthalpy in wellbore
                phases = model.ReservoirModel.getPhaseNames();
                for i = 1:numel(phases)
                    ix = model.ReservoirModel.getPhaseIndex(phases(i));
                    % hWell = feval(model.ReservoirModel.fluid.(['h', phases(i)]), cp, Tw);
                    % h{ix}(injector) = hWell(injector);
                    uWell = feval(model.ReservoirModel.fluid.(['u', phases(i)]), cp, Tw);
                    u{ix}(injector) = uWell(injector);
                end
            end
            % Compute advective heat flux in/out of wellbore
            % q = cellfun(@(rho,h,v) rho.*v.*h, rho, h, v, 'UniformOutput', false);
            q = cellfun(@(rho,u,v) rho.*v.*u, rho, u, v, 'UniformOutput', false);
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