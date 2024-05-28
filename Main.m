function [d, Epsilon, K_hat, Labels, Th, Time, NumSimplices, K1s, WLAPD, BLAPD,PercentKept,Cutoff, knnDistances, k2] = Main(X, Options)

arguments
    X 
    Options.d
    Options.K
    Options.KNNBandWidth
    Options.Noise
    Options.NumScales
    Options.Parallel
    Options.ClusteringMethod
end
if isfield(Options, 'd'), d = Options.d; end
if isfield(Options, 'NumClusters'), K = Options.K; end
if isfield(Options, 'KNNBandWidth'), KNNBandWidth = Options.KNNBandWidth; end
if isfield(Options, 'Noise'), Noise = Options.Noise; end
if isfield(Options, 'NumScales'), NumScales = Options.NumScales; end
if isfield(Options, 'Parallel'), Parallel = Options.Parallel; end
if isfield(Options, 'ClusteringMethod'), ClusteringMethod = Options.ClusteringMethod; end


addpath(genpath(pwd))

if ~exist('Parallel', 'var'), Parallel = 0; end
%if Parallel, parpool; end

if ~exist('d','var') 
    tX = X';     
    EstDimOpts = struct('NumberOfTrials',3,'verbose',0,'MAXDIM',size(X,2),'MAXAMBDIM',size(X,2),'Ptwise',0,'PtIdxs',0,'NetsOpts',[],'UseSmoothedS',0,'EnlargeScales',0,'Deltas',[],'KeepV',0,'kNN',[],'AllPtsFast',0,'DownSample',1,'RandomizedSVDThres',inf);
    fprintf('Intrinsic dimension is not provided. Now estimating... \n');
    [d,Noise_Est] = EstDim_MSVD(tX, Parallel, EstDimOpts);
    fprintf('The intrinsic dimension and noise level found are resp. %d and %.5f. \n', d, Noise_Est);
    if ~exist('Noise','var'), Noise = Noise_Est; end
    clear tX
end

if ~exist('Noise','var') 
    tX = X';     
    EstDimOpts = struct('NumberOfTrials',3,'verbose',0,'MAXDIM',size(X,2),'MAXAMBDIM',size(X,2),'Ptwise',0,'PtIdxs',0,'NetsOpts',[],'UseSmoothedS',0,'EnlargeScales',0,'Deltas',[],'KeepV',0,'kNN',[],'AllPtsFast',0,'DownSample',1,'RandomizedSVDThres',inf);
    fprintf('Noise level is not provided. Now estimating... \n');
    [~,Noise_Est] = EstDim_MSVD(tX, Parallel, EstDimOpts);
    fprintf('The noise level found is %.5f. \n', Noise_Est);
    Noise = Noise_Est;
    clear tX
end

if ~exist('SpectralOpts','var')
    SpectralOpts = struct('NumEig',20,'NumReplicates',5,'RowNormalization',0,'SigmaScaling','Automatic','Laplacian','Symmetric');   
end
if ~isfield(SpectralOpts, 'NumEig'), SpectralOpts.NumEig = 20; end
if ~isfield(SpectralOpts, 'NumReplicates'), SpectralOpts.NumReplicates = 5; end
if ~isfield(SpectralOpts, 'RowNormalization'), SpectralOpts.RowNormalization = 0; end
if ~isfield(SpectralOpts, 'SigmaScaling'), SpectralOpts.SigmaScaling = 'Automatic'; end
if ~isfield(SpectralOpts, 'Laplacian'), SpectralOpts.Laplacian = 'Symmetric'; end


if ~exist('KNNBandWidth','var')
    if d==1, KNNBandWidth=50; C=5;
    elseif d==2, KNNBandWidth = 18; C=7.5; 
    elseif d==3, KNNBandWidth = 38; C=7.5;  
    elseif d==4, KNNBandWidth = 15; C=6;
    end
end  %23, 2.4, 6 %18, 3.5, 5.5 
if ~exist('NumScales','var'), NumScales=30; end

Start = cputime;
%% Randomly sample landmark data points for large d (d >= 3). 
%if d >= 3
%    N = size(X, 1);  %[IDXs,~] = knnsearch(X,X,'K',100);  
%    SampleIdx = sort(randsample(1:N, 2000)); X1 = X(SampleIdx, :); 
%    NoSampleIdx = setdiff(1:N, SampleIdx); 
    %IDXs2 = IDXs(NoSampleIdx, :); %IDXs1 is the NN for sampled data, IDXs2 is the NN for unsampled data. 
%else
    %X1 = X; 
%end

%% Setting Epsilon, the lower bound on nearest neighbors. 
[n, D] = size(X);  %Epsilon = 0.420; %0.45; %%% Need to tune this!!!!
tau = sqrt(D-d)*Noise; Epsilon = 0.50; %C*tau;  %6*tau;  %0.2883; %2.5*Noise; %6*Noise;   %8.6
[IDXs,Dists] = knnsearch(X,X,'K',n); 

[~, K1s]=max(Dists > Epsilon, [], 2); 
K2s = K1s +KNNBandWidth - 1;
NewIDXs = zeros(n, KNNBandWidth); NewDists = zeros(n, KNNBandWidth); 
for k=1:n
    NewIDXs(k,:) = IDXs(k, K1s(k):K2s(k)); NewDists(k,:) = Dists(k, K1s(k):K2s(k));
end
clear Dists

%% Find KNNBandWidth annulus neighbors of each node. 
%{
NewIDXs = zeros(n, KNNBandWidth); %NewDists = zeros(n, KNNBandWidth); 
for k = 1:n
    K1 = min(find(Dists(k,:) > Epsilon)); 
    K1s(k) = K1; 
    %K2 = K1 +KNNBandWidth - 1; NewIDXs(k,:) = IDXs(k, K1:K2);
    %K1s(k) = K1; NewDists(k,:) = Dists(k, K1:K2);
end
clear Dists
%}

%% Form simplex and calculate LAPD. 
%[Simplices, SharedFaces, PosInSim] = BuildSimplices(n,d,NewIDXs,Parallel);
%[Simplices, SharedFaces, PosInSim] = Mutual1(n,d,NewIDXs,Parallel);
[Simplices,SharedFaces,PosInSim,PercentKept] = Filter2(X,n,d,NewIDXs,NewDists,Parallel);
%[Simplice1,SharedFaces1,PosInSim1,PercentKept1] = Filter3(X,n,d,NewIDXs,NewDists,Parallel);
NumSimplices = size(Simplices,1);
[I,J,W] = Adjacency(X,d,Simplices,SharedFaces,PosInSim,Parallel);
clear NewIDXs SharedFaces PosInSim

%% Check connectedness. 
Gknn = sparse(I, J, W, NumSimplices, NumSimplices); Gknn = max(Gknn, Gknn');
[bins,binsizes] = conncomp(graph(Gknn));
if length(binsizes) > 1
    IdxToKeep = bins == mode(bins); Simplices = Simplices(IdxToKeep,:); 
    Gknn = Gknn(IdxToKeep, IdxToKeep); 
    clear IdxToKeep 
    [I,J,W] = find(Gknn);
end
clear Gknn bins binsizes

%% Create CCmatrix. 
[SortedCCmatrix, Th] = ConnectedComponents(I,J,W,Simplices,NumScales);
clear Simplices 
%k2 = ceil(2*(d+1)/log(2)*log(n)*2); %%% Need to tune this!!!!
%[DenoisedCCmatrix, NewTh, Cutoff, knnDistances, k2] = Denoise(d,n,I,J,W,SortedCCmatrix,Th,NumScales);
[DenoisedCCmatrix, NewTh, Cutoff, knnDistances, k2] = Denoise1(d,n,I,J,W,tau,SortedCCmatrix,Th,NumScales);
clear I J W SortedCCmatrix

%if exist('K', 'var') % if the number of components is given by the user. 
%    K_hat = K; 
%    [Labels,WLAPD,BLAPD] = Label_Dend1(n,d,K,DenoisedCCmatrix,NewTh,IDXs);
%end

if strcmp(ClusteringMethod, "Dendrogram")
    if exist('K', 'var') % if the number of components is given by the user. 
        [K_hat,Labels,WLAPD,BLAPD] = Label_Dend(n,d,DenoisedCCmatrix,NewTh,IDXs,K);
    else
        [K_hat,Labels,WLAPD,BLAPD] = Label_Dend(n,d,DenoisedCCmatrix,NewTh,IDXs); 
    end
elseif strcmp(ClusteringMethod, "Spectral")
    %ReasonableColumns = max(DenoisedCCmatrix) < 100 & max(DenoisedCCmatrix) > 1; 
    %ReasonableCCmatrix = min(find(max(DenoisedCCmatrix) < 100));
    %ReasonableCCmatrix = DenoisedCCmatrix(:,[find(ReasonableColumns), NumScales+1:end]);
    [EigVals,EigVecs,Sigmas] = FastEigensolver(d,NumScales,DenoisedCCmatrix,NewTh,SpectralOpts);
    [K_hat1,Labels1,WLAPD1,BLAPD1] = Label_Spec(n, d, NewTh, EigVals, EigVecs, Sigmas, DenoisedCCmatrix, IDXs, SpectralOpts);
end

Time = cputime-Start; 
if Parallel, delete(gcp); end

end