DATAopts.Shape = 'Two Cuboids';    
DATAopts.Number = [2500, 2500];  DATAopts.AmbDim = 4; 
DATAopts.Angles = [0, pi/6];     DATAopts.NoiseSigma = 0.035;

Accuracies = zeros(10, 1); 
maxWLAPDs = zeros(10, 1); modeWLAPDs = zeros(10, 2); 
minBLAPDs = zeros(10, 1); modeBLAPDs = zeros(10, 1); 
Times = zeros(10, 1); 

d=3;Noise=0.035;
for i=1:10  
    i
    [X, LabelsGT] = simdata(DATAopts, i);
    [IntrinsicDim, Epsilon, K_hat, Labels, Th, Time, NumSimplices,K1s,WLAPD, BLAPD,PercentKept,Cutoff,knnDistances,k2] = Main(X, "d",d,"Noise",Noise,"Parallel",0, "ClusteringMethod","Dendrogram");
    %[IntrinsicDim, Epsilon, K_hat, Labels, Th, Time, NumSimplices,K1s,WLAPD, BLAPD,PercentKept] = Main(X,"Parallel",0, "ClusteringMethod","Dendrogram");
    [OA]= GetAccuracies(Labels, LabelsGT, K_hat);
    OA
    Accuracies(i) = OA; Times(i)=Time;
    WithinLAPDs(i) = WLAPD; %BtwLAPDs(i) = BLAPD;
    Cutoffs(i) = Cutoff;
end
[mean(Accuracies) std(Accuracies) mean(Times)]