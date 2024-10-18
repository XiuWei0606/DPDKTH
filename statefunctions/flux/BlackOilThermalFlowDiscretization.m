classdef BlackOilThermalFlowDiscretization < FlowDiscretization
%Discretization and state function grouping for geothermal flow
    
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function props = BlackOilThermalFlowDiscretization(model)
            % Inherit most of the state functions from FLuxDiscretization
            props = props@FlowDiscretization(model);

            % PreDarcy Correction Factor
            props = props.setStateFunction('FacePreDarcyFactor', FacePreDarcyFactor(model));
            % PreDarcy ComponentPhaseFlux
            props = props.setStateFunction('ComponentPhaseFlux', ComponentPreDarcyPhaseFlux(model));
            % PreDarcy PhaseFlux
            props = props.setStateFunction('PhaseFlux', PreDarcyPhaseFlux(model));

            % Fluid flow transmissibility
            % Thermal fracture permeability (NF and EDFM)
            % props = props.setStateFunction('ThermalFactorNf' , ThermalFactorNf(model));
            % props = props.setStateFunction('ThermalFactorF', ThermalFactorF(model));
            % props = props.setStateFunction('Transmissibility', TemTransmissibility(model));

            % Add rock and fluid (water, oil, gas) thermal conductivities
            props = props.setStateFunction('RockThermalConductivity'           , RockThermalConductivity(model) );
            props = props.setStateFunction('MultiphaseFluidThermalConductivity', MultiphaseFluidThermalConductivity(model));
            
            % Rock heat transmissibility
            % Compute static heat transmissibility
            props = props.setStateFunction('RockHeatTransmissibility', MultiphaseHeatTransmissibility(model, 'hr'));

            % Fluid heat transmissibility
            % Compute static heat transmissibility
            props = props.setStateFunction('FluidHeatTransmissibility', MultiphaseHeatTransmissibility(model, 'hf'));

            % Set conductive, advective, and total heat fluxes
            props = props.setStateFunction('ConductiveHeatFlux', MultiphaseConductiveHeatFlux(model));
            props = props.setStateFunction('AdvectiveHeatFlux' , MultiphaseAdvectiveHeatFlux(model));
            props = props.setStateFunction('HeatFlux'          , HeatFlux(model));
            
            % Set molecular diffusion flux % may not need this
            % props = props.setStateFunction('MolecularDiffusivity'       , MolecularDiffusivity(model));
            % props = props.setStateFunction('MolecularTransmissibility'  , DynamicTransmissibility(model, 'MolecularDiffusivity'));
            % props = props.setStateFunction('ComponentTotalDiffusiveFlux', ComponentTotalDiffusiveFlux(model));
        end
        
        %-----------------------------------------------------------------%
        function [acc, flux, names, types] = componentConservationEquations(fd, model, state, state0, dt)
            [acc, flux, names, types] = componentConservationEquations@FlowDiscretization(fd, model, state, state0, dt);
            % Conductive and advective heat flux
            % flowState = fd.buildFlowState(model, state, state0, dt);
            % diffFlux = model.getProp(flowState, 'ComponentTotalDiffusiveFlux');
            % for i = 1:numel(flux)
            %     flux{i} = flux{i} + diffFlux{i};
            % end
        end
        
        %-----------------------------------------------------------------%
        function [acc, flux, name, type] = energyConservationEquation(fd, model, state, state0, dt)
            % Thermal energy in the fluid
            energy = model.getProps(state, 'TotalThermalEnergy');
            % Thermal energy in the fluid at previous timestep
            energy0 = model.getProps(state0, 'TotalThermalEnergy');
            % Conductive and advective heat flux
            flowState = fd.buildFlowState(model, state, state0, dt);
            flux      = model.getProp(flowState, 'HeatFlux');
            % Add up accumulation
            acc = (energy  - energy0 )./dt;
            % Convert to cell arrays to comply with AD framework
            acc  = {acc};
            flux = {flux};
            % Name and type
            name = 'energy';
            type = 'cell';
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