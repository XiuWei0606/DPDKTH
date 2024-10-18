function avg = faceValue2cellValue(G, facevalue)
%Transform face value to cell value by average (exclude zero value)
%
% SYNOPSIS:
%   veclocity = faceFlux2cellVelocity(G, faceFlux)
%
% PARAMETERS:
%   G        - Grid structure.
%
%   faceFlux - Vector of fluxes corresponding to face ordering.
%
% RETURNS:
%   velocity - G.cells.num-by-d matrix of cell velocities.
%
% SEE ALSO:
%   `cellFlux2faceFlux`, `faceFlux2cellFlux`.

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


   % Compute constant part of velocity field in each cell
   if isfield(G, 'FracGrid')
       facevalue = facevalue(1 : length(facevalue) - length(G.nnc.T));
       [cellNo, cellFaces] = getCellNoFaces(G.Matrix);
   else
       [cellNo, cellFaces] = getCellNoFaces(G);
   end
   pick = abs(facevalue) > 0;
   % [cellNo, cellFaces] = getCellNoFaces(G);
   % [cellNo, cellFaces] = getCellNoFaces(G);
   pick = pick(cellFaces);
   a = sparse(cellNo(pick), cellFaces(pick), 1, G.cells.num, G.faces.num) * [facevalue, ones([G.faces.num, 1])];
   avg = a(:,1) ./ a(:,2);
   avg(isnan(avg)) = 0;


   % if ~isfield(G.cells, 'centroids')
   %    G = computeGeometry(G);
   % end
   % 
   % [cellNo, cellFaces] = getCellNoFaces(G);
   % 
   % C  = G.faces.centroids(cellFaces, :) - ...
   %      G.cells.centroids(cellNo   , :);
   % % C(:) = 1;
   % 
   % N = getNeighbourship(G, 'Topological', true);
   % sgn = 2*(N(cellFaces(:,1), 1) == cellNo) - 1;
   % cellvalue = bsxfun(@times, sgn, facevalue(cellFaces, :)); % cell value need to rearrangeds
   % 
   % v = sparse(cellNo, 1 : numel(cellNo), cellvalue) * C;

end
