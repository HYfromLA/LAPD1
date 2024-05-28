function [I,J,W] = Adjacency(X, d, Simplices, SharedFaces, PosInSim, Parallel)

N = size(Simplices, 1); nn = size(SharedFaces, 1); 
%I=[]; J=[]; v1=[]; v2=[];

if Parallel
    parfor i = 1:nn
        SharedNodes = SharedFaces(i, :); VecRows = PosInSim{i}; m=length(VecRows);
        Block = Simplices(VecRows, :); Block = Block';
        OtherNodes = setdiff(Block, SharedNodes,'stable'); UniqueNodes = cat(1, SharedNodes', OtherNodes); 
        coords = X(UniqueNodes,:); temp = (coords - coords(1,:))'; temp = temp(:, 2:end);

        [V1, V2] = PairVectors(d, m, temp);

        TempMat = repmat((1:m)', 1, m); Lower = tril(TempMat, -1); 
        Lower = Lower(:); Lower = Lower(Lower>0); cols = VecRows(Lower); 
        Cumsum = linspace(0, (m-1)*m, m); RepCumsum = repelem(Cumsum, (m-1):-1:0)'; 
        SaveRep2 = Lower+RepCumsum; V2 = V2(SaveRep2,:);
        v1{i} = V1; v2{i} = V2;

        rows = repelem(VecRows(1:(m-1)), (m-1):-1:1, 1);
        I{i} = rows; J{i} = cols; 
    end
else
    for i = 1:nn
        SharedNodes = SharedFaces(i, :); VecRows = PosInSim{i}; m=length(VecRows);
        Block = Simplices(VecRows, :); Block = Block';
        OtherNodes = setdiff(Block, SharedNodes,'stable'); UniqueNodes = cat(1, SharedNodes', OtherNodes); 
        coords = X(UniqueNodes,:); temp = (coords - coords(1,:))'; temp = temp(:, 2:end);

        [V1, V2] = PairVectors(d, m, temp);

        TempMat = repmat((1:m)', 1, m); Lower = tril(TempMat, -1); 
        Lower = Lower(:); Lower = Lower(Lower>0); cols = VecRows(Lower); 
        Cumsum = linspace(0, (m-1)*m, m); RepCumsum = repelem(Cumsum, (m-1):-1:0)'; 
        SaveRep2 = Lower+RepCumsum; V2 = V2(SaveRep2,:);
        v1{i} = V1; v2{i} = V2;
        %v1=[v1;V1]; v2=[v2; V2];

        rows = repelem(VecRows(1:(m-1)), (m-1):-1:1, 1);
        I{i} = rows; J{i} = cols; 
        %I=[I;rows]; J=[J;cols];
    end
end

v1 = cat(1, v1{:}); v2 = cat(1, v2{:}); I = cat(1, I{:}); J = cat(1, J{:}); 
Thetas = real(acos(dot(v1, v2, 2))); 
clear v1 v2

%{
W=pi-Thetas; 
W(W < 1e-8) = 1e-8; 
Keep = W <= pi/2; I = I(Keep); J=J(Keep); W=W(Keep);
%}

W = min(pi-Thetas, Thetas);   
W(W < 1e-8) = 1e-8; 

%A = sparse(I, J, W, N, N); A = max(A, A');

end