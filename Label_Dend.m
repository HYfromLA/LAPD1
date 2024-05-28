function [K_hat, Labels_N, WLAPD, BLAPD] = Label_Dend(n,d,CCmatrix,NewTh,IDXs,K_hat)

%% Some preprocessing. 

m = size(CCmatrix, 2)-(d+1); 

%% Denoise the denoisedCCmatrix if some components contain too few elements, 
%  i.e. < 0.01 of the largest component. 
%Last = max(CCmatrix(:,1:m)); 
Denoised_CCNums = []; WhichColumn = [];k=0;
m1 = max(find(max(CCmatrix(:,1:m))>1)); 
for j=ceil(0.5*m1):m1
    All_CCs = CCmatrix(:, j);
    [Cnts, Uniques] = hist(All_CCs, unique(All_CCs));
    ComponentSize = Cnts ./ sum(Cnts);
    %plot(1:length(ComponentSize), sort(ComponentSize),'.')
    GoodComponents = Uniques(ComponentSize >= 0.1); %% ?
    if length(GoodComponents)>1
        k=k+1;
        Denoised_CCNums(k) = length(GoodComponents); 
        Idx = ismember(All_CCs,GoodComponents); Idxs{k} = Idx; 
        Denoised_CCs{k} = All_CCs(Idx); 
        WhichColumn(k) = j;
    end
end

[Cnts, Uniques] = hist(Denoised_CCNums, unique(Denoised_CCNums)); 

if exist('K_hat')  % If K_hat is given, two cases: (1) a matching in Denoised_CCNums (2) no matching in Denoised_CCNums
    Matches = find(Denoised_CCNums == K_hat);
    if ~isempty(Matches)
        Optimal = max(Matches);
    else

    end
else
    K_hat = Uniques(Cnts == max(Cnts));  
    if length(K_hat) > 1, K_hat = min(K_hat); end 
    Optimal = max(find(Denoised_CCNums == K_hat)); 
end

Labels_S= Denoised_CCs{Optimal}; RemainingSimplices = CCmatrix(Idxs{Optimal}, m+1:end); 
Uniques = unique(Labels_S); for i=1:length(Uniques), Labels_S(Labels_S==Uniques(i))=i; end
Uniques= unique(Labels_S, "stable"); New_Labels_S = zeros(size(Labels_S)); 
for i=1:length(Uniques)
    New_Labels_S(Labels_S==Uniques(i))=i;
end
Labels_S = New_Labels_S; 

%% Find Within LAPD and Between LAPD
RemainingCCmatrix = CCmatrix(Idxs{Optimal}, 1:m);
for i=1:length(Uniques) 
    Component = RemainingCCmatrix(Labels_S==Uniques(i),:); 
    for k=1:m
        if length(unique(Component(:,k))) == 1
            WLAPD(i) = NewTh(k+1); 
            break
        end
    end
end
WLAPD = max(WLAPD); 

BLAPD = NewTh(min(find(max(RemainingCCmatrix)==1))+1);

%% Majority Vote for unique nodes contained in DenoisedCCmatrix. 
RemainingNodes = unique(RemainingSimplices); RemovedNodes = setdiff(1:n, RemainingNodes); 
Labels_N = zeros(n, 1); 
for i=1:length(RemainingNodes)
    [ia, ~] = find(RemainingSimplices == RemainingNodes(i));
    Labels_N(RemainingNodes(i)) = mode(Labels_S(ia));
end

%% Label the removed nodes with majority vote from their 5 nearest labeled neighbors. 
for i = 1:length(RemovedNodes)
    NNs = IDXs(RemovedNodes(i), :); LabeledNeighbors = []; j = 0; 
    while length(LabeledNeighbors) < 5 && j<50
        j=j+1; 
        if Labels_N(NNs(j)) > 0
            LabeledNeighbors = [LabeledNeighbors; NNs(j)];
        end
    end
    Labels_N(RemovedNodes(i)) = mode(Labels_N(LabeledNeighbors));
end

end