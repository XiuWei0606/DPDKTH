function s = setupThermalDPDKEDFMOperatorsTPFA(Gf, rockf, Gm, rockm, tol, varargin)
% Modified from setupEDFMOperatorsTPFA.m.

opt = struct('deck', [], 'neighbors', [], 'trans', [], 'porv', [], ...
             'nfSpacing', [], 'fraction', []);
opt = merge_options(opt, varargin{:});

% tansmissbility and connection in natural fracture grid, including:
% natural fracture grid - natural fracture grid
% EDFM fracture grid - EDFM fracture grid
% natural fracture grid - EDFM fracture grid (NNC type Ⅰ)
% EDFM fracture grid - EDFM fracture grid (NNC type Ⅱ)
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
transmultf = transmultEDFM(Gf, tol);
Tf = Tf.*transmultf;
Thrf = Thrf.*transmultf;
Thff = Thff.*transmultf;
Tf = [Tf; Gf.nnc.T];
Thrf = [Thrf; Gf.nnc.Thr];
Thff = [Thff; Gf.nnc.Thf];
Nf = Gf.faces.neighbors;
Nf =[Nf; Gf.nnc.cells];
intInxf = all(Nf ~= 0, 2);
Tf_all = Tf;
Thrf_all = Thrf;
Thff_all = Thff;
Nf = Nf(intInxf, :); % connection relationship
Tf = Tf(intInxf); % corresponding flowing transmissbility
Thrf = Thrf(intInxf); % corresponding rock heat transmissbility
Thff = Thff(intInxf, :); % corresponding fluid heat transmissbility

% tansmissbility and connection of matrix grid, including:
% matrix grid - matrix grid
% matrix grid - EDFM fracture grid (NNC type Ⅰ)
Tm = getFaceTransmissibility(Gm.Matrix, Gm.Matrix.rock);
lambdaRm = Gm.Matrix.rock.lambdaR.*(1 - Gm.Matrix.rock.poro);
Thrm = getFaceTransmissibility(Gm.Matrix, struct('perm', lambdaRm));
nph = size(Gm.Matrix.rock.lambdaF, 2);
Thfm = zeros(length(Tm), nph);
for i = 1 : nph
    lambdaFm = Gm.Matrix.rock.lambdaF(:, i).*Gm.Matrix.rock.poro;
    Thfm_ = getFaceTransmissibility(Gm.Matrix, struct('perm', lambdaFm));
    Thfm(:, i) = Thfm_;
end
Nm = Gm.Matrix.faces.neighbors + Gf.cells.num;
intInxm = all(Nm ~= Gf.cells.num, 2);
Tm_all = Tm;
Thrm_all = Thrm;
Thfm_all = Thfm;
Nm = Nm(intInxm, :);
Tm = Tm(intInxm);
Thrm = Thrm(intInxm);
Thfm = Thfm(intInxm, :);
% Findex = strcmp(Gf.nnc.type, 'fracmat interior');
% Nm =[Nm; [Gm.nnc.cells(:, 1) + Gf.cells.num Gf.nnc.cells(Findex, 2)]]; % connection relationship
Nm =[Nm; [Gm.nnc.cells(:, 1) + Gf.cells.num Gm.nnc.cells(:, 2)]]; % connection relationship
Tm = [Tm; Gm.nnc.T]; % corresponding flowing transmissbility
Thrm = [Thrm; Gm.nnc.Thr]; % corresponding flowing transmissbility
Thfm = [Thfm; Gm.nnc.Thf]; % corresponding flowing transmissbility
Tm_all = [Tm_all; Gm.nnc.T];
Thrm_all = [Thrm_all; Gm.nnc.Thr];
Thfm_all = [Thfm_all; Gm.nnc.Thf];

% Tm = getFaceTransmissibility(Gm, rockm);
% lambdaRm = rockm.lambdaR.*(1 - rockm.poro);
% Thrm = getFaceTransmissibility(Gm, struct('perm', lambdaRm));
% nph = size(rockm.lambdaF, 2);
% Thfm = zeros(length(Tm), nph);
% for i = 1 : nph
%     lambdaFm = rockm.lambdaF(:, i).*rockm.poro;
%     Thfm_ = getFaceTransmissibility(Gm, struct('perm', lambdaFm));
%     Thfm(:, i) = Thfm_;
% end
% Nm = Gm.faces.neighbors + Gf.cells.num;
% intInxm = all(Nm ~= Gf.cells.num, 2);
% Tm_all = Tm;
% Thrm_all = Thrm;
% Thfm_all = Thfm;
% Nm = Nm(intInxm, :);
% Tm = Tm(intInxm);
% Thrm = Thrm(intInxm);
% Thfm = Thfm(intInxm, :);


% transmissibility between natural fracture grid and matrix grid
Nfm = [(1 : Gf.Matrix.cells.num)' (1 : Gm.Matrix.cells.num)' + Gf.cells.num]; % connection relationship
L = opt.nfSpacing; % natural fracture spacing
fraction = opt.fraction; % volume fractions of natural fracture and matrix
y = @(x) proximity_function(x, L) - fraction(1);
x0 = 0;
xf = fzero(y, x0);
df = xf/2;
Afm = Gm.Matrix.cells.volumes.*minc_proximity_function_dev(xf, L);
dm = getInnerDistance(xf, L);
Tff = Afm.*Gf.Matrix.rock.perm./df;
Tmm = Afm.*Gm.Matrix.rock.perm./dm;
Tfm = 1./(1./Tff + 1./Tmm);
Thrff = Afm.*(1 - Gf.Matrix.rock.poro).*Gf.Matrix.rock.lambdaR./df;
Thrmm = Afm.*(1 - Gm.Matrix.rock.poro).*Gm.Matrix.rock.lambdaR./dm;
Thrfm = 1./(1./Thrff + 1./Thrmm);
Thfff = Afm.*Gf.Matrix.rock.poro.*Gf.Matrix.rock.lambdaF./df;
Thfmm = Afm.*Gm.Matrix.rock.poro.*Gm.Matrix.rock.lambdaF./dm;
Thffm = 1./(1./Thfff + 1./Thfmm);

% assemble all transmissbility and connection relationship
T = [Tf; Tm; Tfm]; % total corresponding flowing transmissbility
Thr = [Thrf; Thrm; Thrfm]; % total corresponding rock heat transmissbility
Thf = [Thff; Thfm; Thffm]; % total corresponding fluid heat transmissbility
N = [Nf; Nm; Nfm]; % total connection relationship
T_all = [Tf_all; Tm_all; Tfm];
Thr_all = [Thrf_all; Thrm_all; Thrfm];
Thf_all = [Thff_all; Thfm_all; Thffm];
numNf = size(Nf, 1);
numNm = size(Nm, 1);
numNfm = size(Nfm, 1);

% poro volume
pvf = poreVolume(Gf, rockf);
pvm = poreVolume(Gm.Matrix, Gm.Matrix.rock);
pv = [pvf; pvm];

% C
nc = Gf.cells.num + Gm.Matrix.cells.num;
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
s.numNm = numNm;
s.numNfm = numNfm;
s.pv = pv;
s.C = C;
s.Grad = @(x) -C*x;
s.Div  = @(x) C'*x;
s.AccDiv = @(acc, flux) acc + C'*flux;
s.faceAvg = @(x) M*x;
s.faceUpstr = upw;
s.splitFaceCellValue = c2f;
s.nfSpacing = opt.nfSpacing;
s.fraction = opt.fraction;
s.Afm = Afm;
s.dm = dm;

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
function [p, dp] = proximity_function(x, L)
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
            dp = 2/L(1);
        case 2
            u = 2*x/L(1);
            v = 2*x/L(2);
            if x < min_x
                p = u + v - u*v;
            else
                p = 1;
            end
            dp = (-8*x + 2*L(1) + 2*L(2))/(L(1)*L(2));
        case 3
            u = 2*x/L(1);
            v = 2*x/L(2);
            w = 2*x/L(3);
            if x < min_x
                p = u + v + w - u*v - u*w - v*w + u*v*w;
            else
                p = 1;
            end
            dp = (24*x^2 - 8*(L(1)+L(2)+L(3))*x +...
                2*(L(1)*L(2) + L(1)*L(3) + L(2)*L(3)))...
                /(L(1)*L(2)*L(3));
        otherwise
            error('Too many entries for fracture spacing L');
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


function D = getInnerDistance(x, L)
    dimension = length(L);
    switch dimension
        case 1
            u = L(1) - 2*x;
            l = u;
            D = l/6;
            % D = 2*l/pi^2;
        case 2
            u = L(1) - 2*x;
            v = L(2) - 2*x;
            l = 2*u*v/(u + v);
            D = l/8;
        case 3
            u = L(1) - 2*x;
            v = L(2) - 2*x;
            w = L(3) - 2*x;
            l = 3*u*v*w/(u*v + v*w + u*w);
            D = l/10;
    end
end