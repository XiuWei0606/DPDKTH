classdef ComponentPreDarcyPhaseFlux < StateFunction
    % Flux of each component, in each phase
    properties (Access = protected)
        mobility_name; % Name of state function where component mobility comes from
    end

    methods
        function cf = ComponentPreDarcyPhaseFlux(model, mob_name)
            cf@StateFunction(model);
            if nargin < 2
                mob_name = 'FaceComponentMobility'; 
            end
            cf.mobility_name = mob_name;
            cf = cf.dependsOn({'PermeabilityPotentialGradient', cf.mobility_name});
            cf = cf.dependsOn('FacePreDarcyFactor');
            cf.label = 'V_{i,\alpha}';
        end

        % Origin
        % function v = evaluateOnDomain(prop, model, state)
        %     ncomp = model.getNumberOfComponents;
        %     nph = model.getNumberOfPhases;
        %     [kgrad, compMob] = prop.getEvaluatedDependencies(state,...
        %         'PermeabilityPotentialGradient', prop.mobility_name);
        %     v = cell(ncomp, nph);
        %     for c = 1:ncomp
        %         for ph = 1:nph
        %             mob = compMob{c, ph};
        %             if ~isempty(mob)
        %                 v{c, ph} = -mob.*kgrad{ph};
        %             end
        %         end
        %     end
        % end

        % PreDarcy
        function v = evaluateOnDomain(prop, model, state)
            ncomp = model.getNumberOfComponents;
            nph = model.getNumberOfPhases;
            [kgrad, compMob] = prop.getEvaluatedDependencies(state,...
                'PermeabilityPotentialGradient', prop.mobility_name);
            v = cell(ncomp, nph);
            if isfield(model.G.rock, 'Mechanisms')
                Coe = prop.getEvaluatedDependencies(state, 'FacePreDarcyFactor');
                for c = 1:ncomp
                    for ph = 1:nph
                        mob = compMob{c, ph};
                        if ~isempty(mob)
                            v{c, ph} = -mob.*kgrad{ph}.*Coe{ph};
                        end
                    end
                end
            else
                for c = 1:ncomp
                    for ph = 1:nph
                        mob = compMob{c, ph};
                        if ~isempty(mob)
                            v{c, ph} = -mob.*kgrad{ph};
                        end
                    end
                end
            end
            % if isfield(model.G.rock, 'THM') && model.G.rock.THM == 1
            %     op = model.operators;
            %     % f
            %     CoeNF = prop.getEvaluatedDependencies(state, 'ThermalFactorNf');
            %     CoeNFface = op.faceAvg(CoeNF);
            %     % F
            %     CoeF = prop.getEvaluatedDependencies(state, 'ThermalFactorF');
            %     CoeFface = op.faceAvg(CoeF);
            % end
        end

        % PreDarcy
        % function v = evaluateOnDomain(prop, model, state)
        %     ncomp = model.getNumberOfComponents;
        %     nph = model.getNumberOfPhases;
        %     [kgrad, compMob] = prop.getEvaluatedDependencies(state,...
        %         'PermeabilityPotentialGradient', prop.mobility_name);
        %     v = cell(ncomp, nph);
        %     if isfield(model.G.rock, 'Mechanisms')
        %         Coe = prop.getEvaluatedDependencies(state, 'FacePreDarcyFactor');
        %         for c = 1:ncomp
        %             for ph = 1:nph
        %                 mob = compMob{c, ph};
        %                 if ~isempty(mob)
        %                     v{c, ph} = -mob.*kgrad{ph}.*Coe{ph};
        %                 end
        %             end
        %         end
        %     else
        %         for c = 1:ncomp
        %             for ph = 1:nph
        %                 mob = compMob{c, ph};
        %                 if ~isempty(mob)
        %                     v{c, ph} = -mob.*kgrad{ph};
        %                 end
        %             end
        %         end
        %     end
        % end


            % 
            % if isfield(model.G.rock, 'Mechanisms')
            %     dp = prop.getEvaluatedDependencies(state, 'PressureGradient');
            %     if isfield(model.G.rock.Mechanisms, 'StaticPreDarcy')
            %         Coe = model.fluid.StaticPreDarcy(dp); % StaticPreDarcy
            %     elseif isfield(model.G.rock.Mechanisms, 'DynamicPreDarcy')
            %         mu = model.getProps(state, 'Viscosity');
            %         for i = 1 : nph
            %             mu{i} = model.operators.faceAvg(mu{i});
            %         end
            %         Coe = model.fluid.DynamicPreDarcy(dp, mu); % StaticPreDarcy
            %     end
            %     if isfield(model.G, 'FracGrid') % PreDarcy only in matrix
            %         a = double(model.G.Matrix.faces.neighbors);
            %         index = all(a ~= 0, 2);
            %         b = a(index, :);
            %         c = [zeros(length(b(:, 1)), 1); ones(length(model.operators.T) - length(b(:, 1)), 1)];
            %         c = logical(c);
            %         for i = 1 : nph
            %             Coe{i}(c) = 1;
            %         end
            %     end
            %     for c = 1:ncomp
            %         for ph = 1:nph
            %             mob = compMob{c, ph};
            %             if ~isempty(mob)
            %                 v{c, ph} = -mob.*kgrad{ph}.*Coe{ph};
            %             end
            %         end
            %     end
            % else
            %     for c = 1:ncomp
            %         for ph = 1:nph
            %             mob = compMob{c, ph};
            %             if ~isempty(mob)
            %                 v{c, ph} = -mob.*kgrad{ph};
            %             end
            %         end
            %     end
            % end


        % PreDarcy
        % function v = evaluateOnDomain(prop, model, state)
        %     ncomp = model.getNumberOfComponents;
        %     nph = model.getNumberOfPhases;
        %     [kgrad, compMob] = prop.getEvaluatedDependencies(state,...
        %         'PermeabilityPotentialGradient', prop.mobility_name);
        %     v = cell(ncomp, nph);
        %     dp = prop.getEvaluatedDependencies(state, 'PressureGradient');
        %     mu = model.getProps(state, 'Viscosity');
        %     hasFracGrid = isfield(G, 'FracGrid');
        %     for i = 1 : nph
        %         mu{i} = model.operators.faceAvg(mu{i});
        %     end
        % 
        %     if isfield(model.G.rock, 'Mechanisms') && isfield(model.G.rock.Mechanisms, 'Nonlinear')
        %         if isfield(model.fluid, 'OilNonlinear')
        %             Coe = model.fluid.OilNonlinear(dp); % 姜
        %         elseif isfield(model.fluid, 'OilNonlinearTemperature')
        %             Coe = model.fluid.OilNonlinearTemperature(dp, mu); % 新
        %         end
        %         % 处理裂缝
        %         a = double(model.G.Matrix.faces.neighbors);
        %         index = all(a ~= 0, 2);
        %         b = a(index, :);
        %         c = [zeros(length(b(:, 1)), 1); ones(length(model.operators.T) - length(b(:, 1)), 1)];
        %         c = logical(c);
        %         % Coe{:}(c) = 1;
        %         for i = 1 : nph
        %              Coe{i}(c) = 1;
        %         end
        %         for c = 1:ncomp
        %             for ph = 1:nph
        %                 mob = compMob{c, ph};
        %                 if ~isempty(mob)
        %                     v{c, ph} = -mob.*kgrad{ph}.*Coe{ph};
        %                 end
        %             end
        %         end
        %     else
        %         for c = 1:ncomp
        %             for ph = 1:nph
        %                 mob = compMob{c, ph};
        %                 if ~isempty(mob)
        %                     v{c, ph} = -mob.*kgrad{ph};
        %                 end
        %             end
        %         end
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
