function [Simplices, SharedFaces, PosInSim, PercentKept] = Filter2(X, n, d, IDXs, Dists, Parallel)
rng('default')

BandWidth = size(IDXs, 2); 


if d==2, Quality_Th = 1.80; 
elseif d==3, Quality_Th= 1.50; %sqrt(2); 
elseif d==4, Quality_Th= 1.80;
end%sqrt(3); %2.1

%IDXs = NewIDXs; Dists = NewDists; 

if d == 1
    IDXs = IDXs'; Simplices = IDXs(:); SharedFaces = (1:n)'; 
    Simplices = cat(2, repelem(SharedFaces, BandWidth, 1), Simplices);
    [Simplices, ~, ~] = unique(sort(Simplices,2),'rows');  
    PercentKept = 1.0;
else
    NodesPerm1 = nchoosek(1:d, 2); NumPerm = nchoosek(d, 2); % Each edge vector has two nodes. %NodesPerm2 = nchoosek(1:(d+1), d);
    m=nchoosek(BandWidth, d);
    for i = 1:n
        NNs = IDXs(i,:); DDs = Dists(i,:);
        Nodes = nchoosek(NNs, d); Node1EdgesNorms{i} = nchoosek(DDs, d);
        
        % All Edges must be long
        a{i} = Nodes(:, NodesPerm1(:,1)); b{i} = Nodes(:, NodesPerm1(:,2));
        %a = Nodes(:, NodesPerm1(:,1)); a = a(:); 
        %b = Nodes(:, NodesPerm1(:,2)); b = b(:);
        %VectorizedEdges = cat(2, a, b);          % All the edges that does not contain node i. i.e. For [1 22 37], this gives [22 37] ([1 22], [1 37] are long enough from KNN and don't need to check).

        %ThirdEdges{i} = VectorizedEdges;
        Simplices{i} = cat(2, repelem(i, m)', Nodes); 
    end
    Node1EdgesNorms = cat(1, Node1EdgesNorms{:});
    a=cat(1, a{:}); b=cat(1, b{:}); a=a(:); b=b(:); ThirdEdges=cat(2,a,b); 
    %ThirdEdges = cat(1, ThirdEdges{:});
    Simplices = cat(1, Simplices{:});
    ThirdEdgesNorms = vecnorm(X(ThirdEdges(:, 2),:) - X(ThirdEdges(:, 1),:), 2, 2);
    ThirdEdgesNorms = reshape(ThirdEdgesNorms, [], NumPerm); EdgesNorms = cat(2, Node1EdgesNorms, ThirdEdgesNorms); % Each row contains the norms of all the edges. 
    Quality = max(EdgesNorms, [], 2) ./ min(EdgesNorms, [], 2); 
    %SortedQuality = sort(Quality); 
    %plot(1:length(SortedQuality),SortedQuality,'.')
    %hold on 
    %xlabel('Index'); ylabel('Max-min Edge Length Ratio');title('Sorted Simplex Qualities')
    %hold off
    %Prompt = 'What quality filter do you want? \n'; Quality_Th = input(Prompt);
    
    GoodQuality = Quality < Quality_Th;
    Simplices = Simplices(GoodQuality, :); 

    % Keep only unique simplices. 
    PercentKept = sum(GoodQuality) / length(GoodQuality); 
    Simplices = unique(sort(Simplices, 2), 'rows'); 
end
%Simplices = Simplices(randsample(1:size(Simplices,1), 0.5*size(Simplices,1)),:);

NodesPerm2 = nchoosek(1:(d+1), d);
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
%Save = (Groups(:, 2)-Groups(:, 1)+1) >= 5;

MoreThanOneAppearance = Groups(Save, :); 
Shared = Shared(ia,:); SharedFaces = Shared(MoreThanOneAppearance(:,1), :);

for i = 1:size(MoreThanOneAppearance, 1)
    PosInSim{i} = sort(Positions(MoreThanOneAppearance(i, 1):MoreThanOneAppearance(i, 2)));
end

end