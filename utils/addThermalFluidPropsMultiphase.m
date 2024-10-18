function fluid = addThermalFluidPropsMultiphase( fluid, varargin )
% Add thermal properties to an existing fluid structure
%
% SYNOPSIS:
%  fluid = addThermalFluidProps(fluid,'pn1', pv1, ...);
%  fluid = addThermalFluidProps(fluid,'pn1', pv1, 'useEOS', true, 'brine', true, ...);
%
% PARAMETERS:
%   fluid   - fluid structure created with initSimpleADIFluid (or other).
%
%   cp      - fluid heat capacity in J kg-1 K-1. typically 4.2e3. Used for
%             the evaluation of the internal energy and enthalpy.
%
%   lambdaF - fluid heat conductivity in W m-1 K-1. typically 0.6.
%
%   useEOS  - logical. By default the density/viscosity formulation
%              of Spivey is used.
% 
%   brine   - logical. Used for the EOS and p,T,c dependency of properties
%             if salt is present or not.
%
%   dNaCl   - salt molecular diffusivity in m2s-1.
% 
% 
%   rho     - density at reservoir condition. Can be given as an handle 
%             function.
% 
%   useBFactor - logical. Used if we want bW definition factor instead 
%                of rhoW 
% 
% RETURNS:
%   fluid - updated fluid structure containing the following functions
%           and properties.
%             * bX(p,T,c)  - inverse formation volume factor
%             * muX(p,T,c) - viscosity functions (constant)
%             * uW(p,T,c)  - internal energy
%             * hW(p,T,c)  - enthalpy
%             * lambdaF    - fluid conductivity
%             * dNaCl      - Salt molecular diffusivity
%
%
% SEE ALSO:
%
% 'initSimpleADIFluid', 'initSimpleThermalADIFluid'

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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

    Watt = joule/second;
    opt = struct('Cp'      , [1, 1, 1]*joule/(kilogram*Kelvin), ...
                 'lambdaF' , [1, 1, 1].*Watt/(meter*Kelvin), ...
                 'EOS'     , false, ...
                 'rho'     , [1, 1, 1].*kilogram/meter^3, ...
                 'cT'      , [0, 0, 0], ...
                 'TRef'    , (273.15 + 20)*Kelvin, ...
                 'pRef'    , 0);
    opt = merge_options(opt, varargin{:});
    % Get phase names
    fNames = fieldnames(fluid);
    phases = '';
    for fNo = 1 : numel(fNames)
        fn = fNames{fNo};
        if strcmpi(fn(1 : 2), 'mu')
            phases = [phases, fn(3)];
        end
    end
    % Set viscosity and density  % XW
    nPh = numel(phases);
    names = upper(phases);
    for phNo = 1 : nPh 
        n = names(phNo);
        if opt.EOS % Calculation based on fitting equations or other functions
            frho = @(varargin) computeEOSDensity(n, varargin{:});
            fluid.(['rho', n]) = frho;
            fluid.(['rho', n, 'S']) = fluid.(['rho', n])(1*atm, opt.TRef);
            fmu = @(varargin) computeEOSViscosity(n, varargin{:});
            fluid.(['mu', n]) = fmu;
        else % only related to pressure
            frho = @(varargin) computeDensity(fluid, n, varargin{:});
            fluid.(['rho', n]) = frho;
        end
    end
    % Set thermal conductivity
    lambdaF = opt.lambdaF;
    fluid.lambdaF = lambdaF;
    % for phNo = 1 : nPh
    %     n = names(phNo);
    %     fluid.(['lambdaF', n]) = lambdaF(phNo);
    % end
    % Set specific, heat capacity, internal energy and enthalpy
    Cp = opt.Cp;
    for phNo = 1 : nPh
        n = names(phNo);
        fluid.(['Cp', n]) = @(varargin) Cp(phNo);
        fluid.(['u', n])  = @(varargin) computeInternalEnergy(fluid, n, varargin{:});
        fluid.(['h', n])  = @(varargin) computeEnthalpy(fluid, n, varargin{:});
    end
end

% Some different functions
%-------------------------------------------------------------------------%
function frho = computeEOSDensity(n, varargin) % XW
    % one can define own density functions for water, oil, gas with
    % different temperature and pressure
    p = varargin{1};
    T = varargin{2};
    if n == 'W'
        frho = density_pure_water(p, T); % kg/m^3
        % frho = 5077.14 - 54.02.*T + 0.28.*T.^2 - 7.14e-4.*T.^3 + 8.81e-7.*T.^4 - 4.26e-10.*T.^5;
        % frho = 1000 + p - p; % kg/m^3
        % frho = 0.00036.*T.^3 - 0.3693.*T.^2 + 122.*T - 0.333.*(p./1e6).^2 + 32.54.*(p./1e6) - 12720; % kg/m^3
    elseif n == 'O'
        % frho = 808.602363 + 6.733190.*(p./1e6).^0.483569 - 0.394209.*(T - 273.15).^1.167098 ...
        %        + 0.034570.*(p./1e6).^0.483569.*(T - 273.15).^1.167098; % kg/m^3
        % frho = 800 + p - p; % kg/m^3
        frho = 800.*(1 + 13.3e-10.*(p - 1*atm) - 9.5e-4.*(T - 293.15 - 20)); % kg/m^3
        % frho = density_pure_water(p, T); % kg/m^3
    elseif n == 'G'
        % frho = 808.602363 + 6.733190.*(p./1e6).^0.483569 - 0.394209.*(T - 273.15).^1.167098 ...
        %        + 0.034570.*(p./1e6).^0.483569.*(T - 273.15).^1.167098; % kg/m^3
        frho = 0.00036.*T.^3 - 0.3693.*T.^2 + 122.*T - 0.333.*(p./1e6).^2 + 32.54.*(p./1e6) - 12720; % kg/m^3
    end
end

%-------------------------------------------------------------------------%
function fmu = computeEOSViscosity(n, varargin) % XW
    % one can define own viscosity functions for water, oil, gas with
    % different temperature and pressure
    p = varargin{1};
    T = varargin{2};
    if n == 'W'
        fmu = viscosity_pure_water(p, T); % Pa*s
        % fmu = 1.233e-4 + 1.466.*exp(-0.025.*T.^5); % Pa*s
        % fmu = 1.233e-4 + 1.466.*exp(-0.025.*T); % Pa*s
        % fmu = 1.233e-4 + 1.466.*exp(-0.025.*T); % Pa*s
        % fmu = 7.14e-9.*T.^2 - 5.642e-6.*T - 5.71e-9.*(p./1e6).^2 + 2.186e-6.*(p./1e6) + 0.0011; % Pa*s
    elseif n == 'O'
        % fmu = 1e-3.*(15.624723 + 0.038356.*(p./1e6) + 0.000253.*(p./1e6).^2 - 0.005130.*(T - 273.15) - 1.133678e-6.*(T - 273.15).^2)./ ...
        %       (1 + 0.008346.*(p./1e6) + 5.19913e-5.*(p./1e6).^2 + 0.005889.*(T - 273.15) + 2.717215e-5.*(T - 273.15).^2); % Pa*s
        % fmu = 1e-3.*34.781.*exp(-0.015.*(T - 273.15));
        % fmu = 1e-3.*17.781.*exp(-0.025.*(T - 273.15)); % Pa*s
        fmu = 1e-3.*34.781.*exp(-0.019.*(T - 273.15)); % Pa*s
        % fmu = viscosity_pure_water(p, T); % Pa*s
    elseif n == 'G'
        % fmu = 1e-3.*(1.624723 + 0.038356.*(p./1e6) + 0.000253.*(p./1e6).^2 - 0.005130.*(T - 273.15) - 1.133678e-6.*(T - 273.15).^2)./ ...
        %       (1 + 0.008346.*v + 5.19913e-5.*(p./1e6).^2 + 0.005889.*(T - 273.15) + 2.717215e-5.*(T - 273.15).^2); % Pa*s
        fmu = 7.14e-9.*T.^2 - 5.642e-6.*T - 5.71e-9.*(p./1e6).^2 + 2.186e-6.*(p./1e6) + 0.0011; % Pa*s
    end
end

%-------------------------------------------------------------------------%
function frho = computeDensity(fluid, n, varargin)
    p = varargin{1};
    T = varargin{2};
    % frho = fluid.(['rho', n, 'S'])+ 0.*p + 0.*T;
    frho = fluid.(['rho', n, 'S']).*fluid.(['b', n])(p);
end

%-------------------------------------------------------------------------%
function u = computeInternalEnergy(fluid, phase, varargin)
    % Cp = feval(fluid.(['Cp', phase]), varargin{:});
    T  = varargin{2};
    Cp = 12010.1471 - 80.4072.*T + 0.3098.*T.^2 - 5.3818e-4.*T.^3 + 3.6253e-7.*T.^4;
    u  = Cp.*T;
end

%-------------------------------------------------------------------------%
function h = computeEnthalpy(fluid, phase, varargin)
    rho = feval(fluid.(['rho', phase]), varargin{:});
    u   = feval(fluid.(['u', phase]), varargin{:});
    p   = varargin{1};
    h   = u + p./rho;
end
