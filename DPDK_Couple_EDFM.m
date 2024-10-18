%% DPDK+EDFM
% f —— natural fracture grid
% F —— EDFM grid
% m —— matrix grid

%%
clear; close all; clc;

%% load modules
mrstModule add ad-blackoil ad-core ad-props mrst-gui hfm shale geothermal compositional linearsolvers
mrstVerbose on

%% parameters
% volume fractions of natural fracture (first) and matrix (second)
fraction = [0.05, 0.95]; 
% natural fracture spacing
nf_Spacing = [100 100 100]; 

%% grid
physdim  = [400, 400, 100];  
celldim  = [40, 40, 10]; 
% f and m share the same grid, but with different volume fractions
% when calculate pore volume, volume fractions should be used
% use ntg to calculate real pore volume
% Gf for f
Gf = cartGrid(celldim, physdim); 
Gf = computeGeometry(Gf);
Gf.rock = makeRock(Gf, 500*milli*darcy, 0.3, 'ntg', fraction(1));
Gf.rock = addThermalRockPropsMultiphase(Gf.rock, ...
                                        'lambdaR', 2, ...
                                        'rhoR', 2700, ...
                                        'CpR', 1000, ...
                                        'lambdaF', 0.8690);
% Gm for m
Gm = cartGrid(celldim, physdim); 
Gm = computeGeometry(Gm);
Gm.rock = makeRock(Gm, 0.1*milli*darcy, 0.08, 'ntg', fraction(2));
Gm.rock = addThermalRockPropsMultiphase(Gm.rock, ...
                                        'lambdaR', 2, ...
                                        'rhoR', 2700, ...
                                        'CpR', 1000, ...
                                        'lambdaF', 0.8690);

%% EDFM fracatures
fracplanes(1).points = [30 20 0;
                        90 150 0;
                        90 150 100;
                        30 20 100];
fracplanes(2).points = [10 40 0;
                        150 60 0;
                        150 60 100;
                        10 40 100];
fracplanes(3).points = [30 160 0;
                        90 90 0;
                        90 90 100;
                        30 160 100];
fracplanes(4).points = [110 30 0;
                        190 190 0;
                        190 190 100;
                        110 30 100];
fracplanes(5).points = [110 90 0;
                        320 130 0;
                        320 130 100;
                        110 90 100];
fracplanes(6).points = [50 230 0;
                        200 140 0;
                        200 140 100;
                        50 230 100];
fracplanes(7).points = [40 180 0;
                        200 300 0;
                        200 300 100;
                        40 180 100];
fracplanes(8).points = [80 320 0;
                        250 190 0;
                        250 190 100;
                        80 320 100];
fracplanes(9).points = [230 160 0;
                        340 240 0;
                        340 240 100;
                        230 160 100];
fracplanes(10).points = [240 80 0;
                         290 240 0;
                         290 240 100;
                         240 80 100];
fracplanes(11).points = [170 210 0;
                         310 290 0;
                         310 290 100;
                         170 210 100];
fracplanes(12).points = [170 210 0;
                         310 290 0;
                         310 290 100;
                         170 210 100];
fracplanes(13).points = [220 290 0;
                         380 260 0;
                         380 260 100;
                         220 290 100];
fracplanes(14).points = [230 270 0;
                         280 370 0;
                         280 370 100;
                         230 270 100];
fracplanes(15).points = [310 200 0;
                         360 390 0;
                         360 390 100;
                         310 200 100];
fracplanes(16).points = [230 340 0;
                         390 360 0;
                         390 360 100;
                         230 340 100];
for i = 1 : numel(fracplanes)
     fracplanes(i).aperture = 5e-3;
     fracplanes(i).poro = 0.3;
     fracplanes(i).perm = 2e4*milli*darcy;
     fracplanes(i).lambdaR = 2;
     fracplanes(i).lambdaF = 0.8690;
     fracplanes(i).rhoR = 2700;
     fracplanes(i).CpR = 1000;
     fracplanes(i).ntg = 1;
end
figure;
plotfracongrid(Gf, fracplanes);
view(2)

% deal fracplanes
[fracplanesm, fracplanesf] = deal(fracplanes);

%% EDFM PreProcessing
tol = 1e-5;
% obtain f-F transmissbility: NNC type Ⅰ and Ⅱ 
[Gf, fracplanesf] = EDFMshalegridThermal(Gf, fracplanesf, 'Tolerance', tol, 'fracturelist', 1 : numel(fracplanesf));
Gf = fracturematrixShaleNNC3DThermal(Gf, tol);
[Gf, fracplanesf] = fracturefractureShaleNNCs3DThermal(Gf, fracplanesf, tol, 'Verbose', true);
% obtain m-F transmissbility: NNC type Ⅰ
[Gm, fracplanesm] = EDFMshalegridThermal(Gm, fracplanesm, 'Tolerance', tol, 'fracturelist', 1 : numel(fracplanesm));
Gm = fracturematrixShaleNNC3DThermal(Gm, tol);

%% operators
% obtain total transmissbility
operators = setupThermalDPDKEDFMOperatorsTPFA(Gf, Gf.rock, Gm, Gm.rock, tol, ...
                                              'nfSpacing', nf_Spacing, ...
                                              'fraction', fraction);
% index of f, F, m
grid_f = (1 : Gf.Matrix.cells.num)';
grid_F = (Gf.Matrix.cells.num + 1 : Gf.cells.num)';
grid_m = (Gf.cells.num + 1 : Gf.cells.num + Gm.Matrix.cells.num)';

%% fluid
fluid = initSimpleADIFluid('mu', 1e-3, 'rho', 1000, 'phases', 'W');
% same as fluid function in geothermal module
fluid = addThermalFluidPropsMultiphase(fluid, 'Cp', 4.2e3, ...
                                              'lambdaF', 0.8690, ...
                                              'EOS', true);

%% model
gravity reset off
% similar to Geothermalmodel
model =  BlackOilThermalModel(Gf, Gf.rock, fluid, 'water', true, 'oil', false, 'gas', false);
model.operators = operators;
% total grid number
model.G.cells.num = Gf.cells.num + Gm.Matrix.cells.num;
% total grid centriods
model.G.cells.centroids = [Gf.cells.centroids; Gm.Matrix.cells.centroids];
% real cell volume should be corrected by volume fractions (ntg)
model.G.cells.volumes = [Gf.cells.volumes.*Gf.rock.ntg; Gm.Matrix.cells.volumes.*Gm.Matrix.rock.ntg];
model.operators.vol = model.G.cells.volumes;
% physical parameters should have the same length as total grid number
model.rock = struct('perm', [Gf.rock.perm; Gm.Matrix.rock.perm], ...
                    'poro', [Gf.rock.poro; Gm.Matrix.rock.poro], ...
                    'lambdaR', [Gf.rock.lambdaR; Gm.Matrix.rock.lambdaR], ...
                    'rhoR', [Gf.rock.rhoR; Gm.Matrix.rock.rhoR], ...
                    'lambdaF', [Gf.rock.lambdaF; Gm.Matrix.rock.lambdaF], ...
                    'CpR', [Gf.rock.CpR; Gm.Matrix.rock.CpR], ...
                    'uR', @(p, T) rock.CpR.*T, ...
                    'tau', [Gf.rock.tau; Gm.Matrix.rock.tau]);
% should be 0!
model.outputFluxes = 0;
model.OutputStateFunctions = {'ConductiveHeatFlux', 'AdvectiveHeatFlux'};

%% initial state
state0  = initResSol(model.G, 300*barsa, 1);
state0.T = ones(model.G.cells.num, 1).*(273.15 + 200);

%% wells
injf = Index2Num(Gf, 5, 5, 1 : 10);
prof = Index2Num(Gf, 36, 36, 1 : 10);
massrate = 5*kilogram/second;
rate = massrate/fluid.rhoWS;
W = addWell([], Gf.Matrix, Gf.Matrix.rock, injf, 'type', 'rate', 'Val', rate, ...
            'compi', 1, 'name', 'Inj', 'sign', 1);
W = addWell(W, Gf.Matrix, Gf.Matrix.rock, prof, 'type', 'bhp', 'Val', 300*barsa, ...
            'compi', 1, 'name', 'Inj', 'sign', -1);
W = addThermalWellProps(W, Gf, Gf.rock, fluid, 'T', 303.15);
% ignore heat conduction in well
for i = 1 : numel(W)
    W(i).WIth(:) = 0;
end

figure; 
plotGrid(Gf, 'facealpha', 0, 'edgealpha', 0.1); axis equal tight; hold on;
plotWell(Gf, W, 'height', 0.1); hold on;

%% schedule
dt = rampupTimesteps(30*year, 30*day, 10);
time = cumsum(dt)./year;
schedule = simpleSchedule(dt, 'W', W);

%% simulation
[ws, states, report] = simulateScheduleAD(state0, model, schedule);

% Here is the acceleration problem, bad convergence
% model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true);
% lsolver = AMGCL_CPRSolverAD('tolerance', 1e-4);
% lsolver = GMRES_ILUSolverAD('tolerance', 1e-4);
% [ws, states, report] = simulateScheduleAD(state0, model, schedule, ...
%                                           'linearSolver', lsolver);

%% resutls - curves
plotWellSols(ws, cumsum(schedule.step.val))

%% resutls - distributions
figure('Position', [100, 50, 700, 800]);
subplot(3, 1, 1)
plotCellData(Gf.Matrix, states{132}.T(grid_f, 1) - 273.15, 'edgealpha', 0); axis equal tight;
plotGrid(Gf, grid_F, 'facealpha', 0)
colorbar; colormap(jet); view(3)
clim([20 200]);
subplot(3, 1, 2)
plotCellData(Gf.Matrix, states{132}.T(grid_m, 1) - 273.15, 'edgealpha', 0); axis equal tight;
plotGrid(Gf, grid_F, 'facealpha', 0)
colorbar; colormap(jet); view(3)
clim([20 200]);
subplot(3, 1, 3)
plotCellData(Gf, states{132}.T(grid_F, 1) - 273.15, grid_F, 'edgealpha', 0); axis equal tight;
plotGrid(Gf, grid_F, 'facealpha', 0)
colorbar; colormap(jet); view(3)
clim([20 200]);
