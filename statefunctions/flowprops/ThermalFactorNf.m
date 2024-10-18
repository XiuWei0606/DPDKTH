classdef ThermalFactorNf < StateFunction
%State function for variable factor of NF permeability

    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = ThermalFactorNf(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn('temperature', 'state');
            gp = gp.dependsOn('pressure', 'state');
            gp.label = 'C_NF';
        end
        
        %-----------------------------------------------------------------%
        function CNF = evaluateOnDomain(prop, model, state)
            [p, T] = model.getProps(state, 'pressure', 'temperature');
            CNF = model.fluid.FactorNf(p, T); 
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