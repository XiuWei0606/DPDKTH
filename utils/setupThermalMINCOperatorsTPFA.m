function s = setupThermalMINCOperatorsTPFA(Gf, rockf, GmTotal, varargin)
% Modified from setupEDFMOperatorsTPFA.m.

opt = struct('deck', [], 'neighbors', [], 'trans', [], 'porv', [], ...
             'nfSpacing', [], 'fraction', [], 'region', [], 'MatrixSize', []);
opt = merge_options(opt, varargin{:});

% tansmissbility and connection in natural fracture grid, including:
% natural fracture grid - natural fracture grid
Tf = getFaceTransmissibility(Gf, rockf);
lambdaRf = rockf.lambdaR.*(1 - rockf.poro);
Thrf = getFaceTransmissibility(Gf, struct('perm', lambdaRf));
nph = size(rockf.lambdaF, 2);
Thff = zeros(length(Tf), nph);
for i = 1 : nph
    lambdaFf = rockf.lambdaF(:, i).*rockf.poro;
    Thff_ = getFaceTransmissibility(Gf, struct('perm', lambdaFf));
    Thff(:, i) = Thff_;
end
Nf = Gf.faces.neighbors;
intInxf = all(Nf ~= 0, 2);
Tf_all = Tf;
Thrf_all = Thrf;
Thff_all = Thff;
Nf = Nf(intInxf, :); % connection relationship
Tf = Tf(intInxf); % corresponding flow transmissbility
Thrf = Thrf(intInxf); % corresponding rock heat transmissbility
Thff = Thff(intInxf, :); % corresponding fluid heat transmissbility

% MINC parameters
[x, dist, area] = MINCParameters(opt.fraction, opt.nfSpacing, opt.MatrixSize);
d = zeros(length(opt.fraction), 1);
d(1) = x(1)/2;
d(end) = x(end) - x(end - 1);
for i = 2 : length(d) - 1
    d(i) = (x(i) - x(i - 1))/2;
end

% connection between natural fracture grid and first matrix grid
FirstMatrixIndexMap = Gf.cells.num + GmTotal{1}.cells.indexMap;
Nfm = [opt.region FirstMatrixIndexMap]; % connection between natural fracture grid and first matrix grid

% connection between matrix grid
for i = 2 : numel(GmTotal)
    Nmm = [FirstMatrixIndexMap + (i - 2).*GmTotal{i}.cells.num FirstMatrixIndexMap + (i - 1).*GmTotal{i}.cells.num];
    Nfm = [Nfm; Nmm];
end

% transmissbility
Tff = area(1).*rockf.perm(opt.region)./d(1);
Thrff = area(1).*rockf.lambdaR(opt.region).*(1 - rockf.poro(opt.region))./d(1);
Thfff = area(1).*rockf.lambdaF(opt.region, :).*rockf.poro(opt.region)./d(1);
Tmm1 = area(1).*GmTotal{1}.rock.perm./d(2);
Thrmm1 = area(1).*GmTotal{1}.rock.lambdaR.*(1 - GmTotal{1}.rock.poro)./d(1);
Thfmm1 = area(1).*GmTotal{1}.rock.lambdaF.*GmTotal{1}.rock.poro./d(1);
Tfm = 1./(1./Tff + 1./Tmm1); % flow transmissbility between fracture grid and first matrix grid
Thrfm = 1./(1./Thrff + 1./Thrmm1); % rock heat transmissbility between fracture grid and first matrix grid
Thffm = 1./(1./Thfff + 1./Thfmm1); % fluid heat transmissbility between fracture grid and first matrix grid
% Tfm = Tmm1; % flow transmissbility between fracture grid and first matrix grid
% Thrfm = Thrmm1; % rock heat transmissbility between fracture grid and first matrix grid
% Thffm = Thfmm1; % fluid heat transmissbility between fracture grid and first matrix grid
for i = 2 : numel(GmTotal)
    TmmFormer = area(i).*GmTotal{i - 1}.rock.perm./d(i);
    ThrmmFormer = area(i).*GmTotal{i - 1}.rock.lambdaR.*(1 - GmTotal{i - 1}.rock.poro)./d(i);
    ThfmmFormer = area(i).*GmTotal{i - 1}.rock.lambdaF.*GmTotal{i - 1}.rock.poro./d(i);
    TmmLater =  area(i).*GmTotal{i}.rock.perm./d(i + 1);
    ThrmmLater =  area(i).*GmTotal{i}.rock.lambdaR.*(1 - GmTotal{i}.rock.poro)./d(i + 1);
    ThfmmLater =  area(i).*GmTotal{i}.rock.lambdaF.*GmTotal{i}.rock.poro./d(i + 1);
    % Tmm = 1./(1./TmmFormer + 1./TmmLater);
    % Thrmm = 1./(1./ThrmmFormer + 1./ThrmmLater);
    % Thfmm = 1./(1./ThfmmFormer + 1./ThfmmLater);
    Tmm = 1./(1./TmmFormer + 1./TmmLater);
    Thrmm = 1./(1./ThrmmFormer + 1./ThrmmLater);
    Thfmm = 1./(1./ThfmmFormer + 1./ThfmmLater);
    Tfm = [Tfm; Tmm];
    Thrfm = [Thrfm; Thrmm];
    Thffm = [Thffm; Thfmm];
end

% assemble all transmissbility and connection relationship
T = [Tf; Tfm]; % total flow corresponding transmissbility
Thr = [Thrf; Thrfm]; % total flow corresponding transmissbility
Thf = [Thff; Thffm]; % total flow corresponding transmissbility
N = [Nf; Nfm]; % total connection relationship
T_all = [Tf_all; Tfm];
Thr_all = [Thrf_all; Thrfm];
Thf_all = [Thff_all; Thffm];
numNf = size(Nf, 1);
numNfm = size(Nfm, 1);

%
indexMapf = (1 : Gf.cells.num)';
indexMapm = (Gf.cells.num + 1 : Gf.cells.num + numel(GmTotal)*numel(opt.region))';

% poro volume
pvf = poreVolume(Gf, rockf);
pvm = [];
for i = 1 : numel(GmTotal)
    pvm = [pvm; poreVolume(GmTotal{i}, GmTotal{i}.rock)];
end
pv = [pvf; pvm];

% C
nc = Gf.cells.num + numel(GmTotal)*numel(opt.region);
n = size(N, 1);
C  = sparse([(1 : n)'; (1 : n)'], N, ones(n, 1)*[1 -1], n, nc);

% faceAvg
nf = size(N, 1);
M  = sparse((1 : nf)'*[1 1], N, .5*ones(nf, 2), nf, nc);

% faceUpstr
upw = @(flag, x) faceUpstr(flag, x, N, [nf, nc]);

% splitFaceCellValue
c2f = @(operators, flag, x) splitFaceCellValue(operators, flag, x, [nf, nc]);

% include
s.T = T;
s.Thr = Thr;
s.Thf = Thf;
s.T_all = T_all;
s.Thr_all = Thr_all;
s.Thf_all = Thf_all;
s.N = N;
s.numNf = numNf;
s.numNfm = numNfm;
s.indexMapf = indexMapf;
s.indexMapm = indexMapm;
s.region = opt.region;
s.pv = pv;
s.C = C;
s.Grad = @(x) -C*x;
s.Div  = @(x) C'*x;
s.AccDiv = @(acc, flux) acc + C'*flux;
s.faceAvg = @(x) M*x;
s.faceUpstr = upw;
s.splitFaceCellValue = c2f;

end

function xu = faceUpstr(flag, x, N, sz)
    if numel(flag) == 1
        flag = repmat(flag, size(N, 1), 1);
    end
    assert(numel(flag) == size(N, 1) && islogical(flag), ...
        'One logical upstream flag must'' be supplied per interface.');
    upCell       = N(:,2);
    upCell(flag) = N(flag,1);
    if isnumeric(x)
        % x is a simple matrix, we just extract values using the cells
        xu = x(upCell, :);
    else
        % x is likely AD, construct a matrix to achieve upstream weighting
        xu = sparse((1:sz(1))', upCell, 1, sz(1), sz(2))*x;
    end
end

function [x, dist, area] = MINCParameters(fractions, nfSpacing, MatrixSize)
    p = zeros(length(fractions) - 1, 1);
    v = 0;
    for i = 1 : length(p)
        v = v + fractions(i);
        p(i) = v;
    end
    x = zeros(length(p), 1);
    for i = 1 : length(x)
        y = @(x) minc_proximity_function(x, nfSpacing) - p(i);
        x0 = 0;
        x(i) = fzero(y, x0);
    end
    dist = zeros(length(p), 1);
    dist(1) = x(1)/2 + (x(2) - x(1))/2;
    for i = 2 : length(dist) - 1
        dist(i) = (x(i + 1) - x(i))/2 + (x(i) - x(i - 1))/2;
    end
    D = getMostInnerDistance(x(end), nfSpacing);
    dist(end) = D + (x(end) - x(end -1))/2;
    x(end + 1) = x(end) + D;
    area = zeros(length(p), 1);
    for i = 1 : length(area)
        area(i) = prod(MatrixSize)*minc_proximity_function_dev(x(i), nfSpacing);
    end
end

function p = minc_proximity_function(x, L)
    min_x = min(L);
    dimension = length(L);
    switch dimension
        case 1
            u = 2*x/L(1);
            if x < min_x
                p = u;
            else
                p = 1;
            end
        case 2
            u = 2*x/L(1);
            v = 2*x/L(2);
            if x < min_x
                p = u + v - u*v;
            else
                p = 1;
            end
        case 3
            u = 2*x/L(1);
            v = 2*x/L(2);
            w = 2*x/L(3);
            if x < min_x
                p = u + v + w - u*v - u*w - v*w + u*v*w;
            else
                p = 1;
            end
        otherwise
            error('Too many entries for fracture spacing L');
    end
end

function D = getMostInnerDistance(x, L)
    dimension = length(L);
    switch dimension
        case 1
            u = L(1) - 2*x;
            l = u;
            % D = l/6;
            D = 2*l/pi^2;
            % D = mean([u/4]);
        case 2
            u = L(1) - 2*x;
            v = L(2) - 2*x;
            l = 2*u*v/(u + v);
            % D = l/8;
            D = 2*l/pi^2;
            % D = mean([u/4, v/4]);
        case 3
            u = L(1) - 2*x;
            v = L(2) - 2*x;
            w = L(3) - 2*x;
            l = 3*u*v*w/(u*v + v*w + u*w);
            % D = l/10;
            D = 2*l/pi^2;
            % D = mean([u/4, v/4, w/4]);
    end
end

function dp = minc_proximity_function_dev(x, L)
    dimension = length(L);
    switch dimension
        case 1
            dp = 2/L(1);
        case 2
            dp = (-8*x + 2*L(1) + 2*L(2))/(L(1)*L(2));
        case 3
            dp = (24*x^2 - 8*(L(1)+L(2)+L(3))*x +...
                2*(L(1)*L(2) + L(1)*L(3) + L(2)*L(3)))...
                /(L(1)*L(2)*L(3));
        otherwise
            error('Too many entries for fracture spacing L');
    end
end