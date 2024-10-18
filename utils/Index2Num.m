function Num = Index2Num(G, I, varargin)
assert (isnumeric(I) && numel(I) == 1);

%--------------------------------------------------------------------------
% Determine call mode. ----------------------------------------------------
%
if mod(nargin, 2) == 0
   % W = verticalWell(W, G, rock, I, J, K, ...)
   %
   assert (isfield(G, 'cartDims'));
   layerSize = prod(G.cartDims(1:2));
   if G.griddim == 3
      numLayers = G.cartDims(3);
   else
      numLayers = 1;
   end

   J = varargin{1}; varargin = varargin(2 : end);
   assert (isnumeric(J) && numel(J) == 1);

   assert ((1 <= I) && (I <= G.cartDims(1)), 'I-position out of bounds');
   assert ((1 <= J) && (J <= G.cartDims(2)), 'J-position out of bounds');

   I = I + (J - 1)*G.cartDims(1);    clear J
else
   % W = verticalWell(W, G, rock, I, K, ...)
   %
   if isfield(G, 'layerSize') && isfield(G, 'numLayers')
      % G presumably a 'makeLayeredGrid'.
      layerSize = G.layerSize;
      numLayers = G.numLayers;
   elseif isfield(G, 'cartDims')
      % Logically Cartesian grid.
      layerSize = prod(G.cartDims(1:2));
      if G.griddim == 3
         numLayers = G.cartDims(3);
      else
         numLayers = 1;
      end
   else
      error(msgid('GridType:NotLayered'), ...
            'Grid type does not appear to support a layered structure');
   end
end

%--------------------------------------------------------------------------
% Extract layer vector (K). -----------------------------------------------
%
K = varargin{1}; varargin = varargin(2 : end);
if isempty(K)
   % Empty 'K'.  Complete in all layers.
   K = (1 : numLayers) .';
end
assert ((1 <= min(K)) && (max(K) <= numLayers));

%--------------------------------------------------------------------------
% Extract active completions. ---------------------------------------------
%
act                   = zeros([numLayers * layerSize, 1]);
act(G.cells.indexMap) = 1 : numel(G.cells.indexMap);

wc = act(I + (K - 1).*layerSize);
i  = wc > 0;
Num = wc(i);
end

