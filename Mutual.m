function [Simplices, SharedFaces, PosInSim] = Mutual(n, d, IDXs, Parallel)

rng('default')

BandWidth = size(IDXs, 2); 

%% Build simplices with annulus KNN graph on X.

% All valid vectors (vectors in the annulus KNN graph). 

if Parallel
    parfor i=1:n
        Node_i_nbrs = IDXs(i,:)'; [a,~] = find(IDXs==i); % Mutual KNN.
        Node_i_nbrs = union(Node_i_nbrs, a); m=length(Node_i_nbrs);
        Neighbors{i} = Node_i_nbrs;
        ValidVectors{i} = cat(1, cat(2, repelem(i, m)', Node_i_nbrs), cat(2, Node_i_nbrs,repelem(i, m)')); 
    end
    ValidVectors = unique(cat(1,ValidVectors{:}),"rows");

    NodesPerm = nchoosek(1:d, d-1)';
    parfor i = 1:n
        Node_i_nbrs = Neighbors{i};
        Nodes = nchoosek(Node_i_nbrs, d); 
        if d == 2 
            flag= ismember(Nodes, ValidVectors, "rows"); 
        else
            v = Nodes(:, NodesPerm'); Reshaped_v = [];
            for j = 1 : d 
                Reshaped_v = cat(1, Reshaped_v, v(:, (j-1)*2+1:j*2));     
            end       
            flag_v= ismember(Reshaped_v, ValidVectors, "rows"); 
            step = size(v,1); flag = flag_v(1:step);
            for j=1:(d-1)
                flag = and(flag, flag_v(j*step+1:(j+1)*step));
            end
        end

        Nodes = Nodes(flag,:); Nodes = cat(2, repelem(i, sum(flag))', Nodes);
        Simplices{i} = Nodes; 
    end
else
    for i=1:n
        Node_i_nbrs = IDXs(i,:)'; [a,~] = find(IDXs==i); % Mutual KNN.
        Node_i_nbrs = union(Node_i_nbrs, a); m=length(Node_i_nbrs);
        Neighbors{i} = Node_i_nbrs;
        ValidVectors{i} = cat(1, cat(2, repelem(i, m)', Node_i_nbrs), cat(2, Node_i_nbrs,repelem(i, m)')); 
    end
    ValidVectors = unique(cat(1,ValidVectors{:}),"rows");

    NodesPerm = nchoosek(1:d, d-1)';
    for i = 1:n
        Node_i_nbrs = Neighbors{i};
        Nodes = nchoosek(Node_i_nbrs, d); 
        if d == 2
            flag= ismember(Nodes, ValidVectors, "rows"); 
        else
            v = Nodes(:, NodesPerm); Reshaped_v = [];
            for j = 1 : d 
                Reshaped_v = cat(1, Reshaped_v, v(:, (j-1)*2+1:j*2));     
            end       
            flag_v= ismember(Reshaped_v, ValidVectors, "rows"); 
            step = size(v,1); flag = flag_v(1:step);
            for j=1:(d-1)
                flag = and(flag, flag_v(j*step+1:(j+1)*step));
            end
        end

        Nodes = Nodes(flag,:); Nodes = cat(2, repelem(i, sum(flag))', Nodes);
        Simplices{i} = Nodes; 
    end
end 
% Keep only unique simplices. 
Simplices = cat(1, Simplices{:}); Simplices = unique(sort(Simplices, 2), 'rows'); 


N = size(Simplices, 1); NodesPerm2 = nchoosek(1:(d+1), d);
v = Simplices(:, NodesPerm2');  Positions = repmat((1:size(Simplices, 1))', (d+1), 1);
Shared = []; 
for j = 1 : d 
    w = v(:, j:d:end);  w = w(:); Shared = cat(2,Shared, w);
end

[~, ~, T] = unique(Shared, 'rows');
[a, ia] = sort(T); nn = length(a); Positions = Positions(ia);
groupEnd = [find(diff(a)==1); nn];  groupStart = [1; groupEnd(1:end-1)+1];
Groups = cat(2, groupStart, groupEnd); 
Save = Groups(:, 1) ~= Groups(:, 2);  % Keep only shared faces that appear more than once (so that at least an angle can form).

MoreThanOneAppearance = Groups(Save, :); 
Shared = Shared(ia,:); SharedFaces = Shared(MoreThanOneAppearance(:,1), :);

if Parallel
    parfor i = 1:size(MoreThanOneAppearance, 1)
        PosInSim{i} = sort(Positions(MoreThanOneAppearance(i, 1):MoreThanOneAppearance(i, 2)));
    end
else
    for i = 1:size(MoreThanOneAppearance, 1)
        PosInSim{i} = sort(Positions(MoreThanOneAppearance(i, 1):MoreThanOneAppearance(i, 2)));
    end
end

end