function [sortedCCmatrix, th] = ConnectedComponents(I,J,W,Simplices,NumScales)
  
% ConnectedComponents.m
% 
% IN: 
% 
% Gknn: The weight matrix; caculated in ComputeAdjacency.m.
% I, J, PositiveWeights: positions and values of positive entries in Gknn; caculated in ComputeAdjacency.m.
% Simplices: An nn1 by d matrix; calculated in FindSimplices.m. 
% LAPDOpts: Structure of options for computing LAPD. 
%
% OUT:
%
% sortedCCmatrix: a matrix containing the connected components at each scale
% and containing simplices as the final d+1 columns. This matrix is sorted
% from rightmost column to left for efficient path distance querying.
% th: vector of thresholds corresponding to different scales. 

N = size(Simplices, 1);

%% Determine thresholds.
tmin = prctile(W, 1); %Take 1st percentile as smallest scale
tmax = 1.01*max(W); th = exp(linspace(log(tmin),log(tmax),NumScales)); 
th=[0,th];

%% Threshold graph.
for i=2:(NumScales+1)
    Idx = W <= th(i); R=I(Idx); C=J(Idx); E=W(Idx);
    G = sparse(R, C, E, N, N); G = max(G, G');
    [~, CCmatrix(:,i-1)] = AlternateConnComp(G);
end
clear G

%% Create Sorted Matrix of Component Indices
CCmatrix = cat(2, CCmatrix, Simplices); 
sortedCCmatrix = sortrows(CCmatrix, NumScales:-1:1);

end