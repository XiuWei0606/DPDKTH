classdef FacePreDarcyFactor < StateFunction
    % PreDarcy incorrection factor of each phase in grid face
    properties
        
    end

    methods
        function fc = FacePreDarcyFactor(model)
            fc@StateFunction(model);
            fc = fc.dependsOn({'PhasePotentialDifference'});
            fc = fc.dependsOn('Viscosity', 'PVTPropertyFunctions');
            fc.label = 'Chi_\alpha';
        end

        % PreDarcy
        function Coe = evaluateOnDomain(prop, model, state)
            nph = model.getNumberOfPhases;
            if isfield(model.G.rock, 'Mechanisms')
                dp = prop.getEvaluatedDependencies(state, 'PhasePotentialDifference');
                % dl = model.operators.Grad(model.G.cells.centroids);
                % [col_indices, row_indices] = find(dl');
                % dl_ = dl(sub2ind(size(dl), row_indices, col_indices));
                % dl_ = dl_(1 : size(dl));
                if isfield(model.G.rock.Mechanisms, 'StaticPreDarcy')
                    Coe = model.fluid.StaticPreDarcy(dp); % StaticPreDarcy
                elseif isfield(model.G.rock.Mechanisms, 'DynamicPreDarcy')
                    mu = model.getProps(state, 'Viscosity');
                    for i = 1 : nph
                        mu{i} = model.operators.faceAvg(mu{i});
                    end
                    Coe = model.fluid.DynamicPreDarcy(dp, mu); % StaticPreDarcy
                end
                % if isfield(model.G, 'FracGrid') % PreDarcy only in matrix
                % EDFM + M
                %     a = double(model.G.Matrix.faces.neighbors);
                %     index = all(a ~= 0, 2);
                %     b = a(index, :);
                %     c = [zeros(length(b(:, 1)), 1); ...
                %         ones(length(model.operators.T) - length(b(:, 1)) - length(model.G.nnc.T), 1); ...
                %         zeros(length(model.G.nnc.T), 1);];
                %     c = logical(c);
                %     for i = 1 : nph
                %         Coe{i}(c) = 1;
                %     end
                % end
                if isfield(model.operators, 'numNfm') % PreDarcy only in matrix 
                    % EDFM +DPDK
                    a = model.operators.numNf + model.operators.numNm + model.operators.numNfm;
                    index = ones(a, 1);
                    index(model.operators.numNf + 1 : model.operators.numNf + model.operators.numNm) = 0;
                    index = logical(index);
                    for i = 1 : nph
                        Coe{i}(index) = 1;
                    end
                end
            else
                Coe = cell(1, nph);
                for i = 1 : nph
                    Coe{i} = 1;
                end
            end
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
