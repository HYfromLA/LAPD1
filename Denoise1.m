function [DenoisedCCmatrix, NewTh, Cutoff, Ths, k2] = Denoise1(d, n, I, J, W, tau,sortedCCmatrix, th, NumScales)

rng('default')
nn = size(sortedCCmatrix,1); m = size(sortedCCmatrix,2)-(d+1); 

k2 = floor(2*(d+1)/log(2)*log(n)*2);  

knndistances = zeros(nn, 1); 
    for s=1:nn
        iup = s; %will tell us the index of a nearest neighbor when it is found
        idown = s; %will tell us the index of a nearest neighbor when it is found
        counter = 2; %keeps track of how many nearest neighbors you have found

        for j=1:m %defines the approximate path distance
            % Move up until you find a point not in the same CC then move
            % left and repeat
            while (iup>1 && sortedCCmatrix(iup,j)==sortedCCmatrix(iup-1,j) && counter<k2+1)
                iup = iup-1;
                if counter==k2
                    knndistances(s)=th(j+1);
                end
                counter=counter+1;
            end
            % Move down until you find a point not in the same CC then move
            % left and repeat
            while (idown<nn && sortedCCmatrix(idown,j)==sortedCCmatrix(idown+1,j) && counter<k2+1)
                idown = idown+1;
                if counter==k2
                    knndistances(s)=th(j+1);
                end
                counter=counter+1;
            end
        end
        % If we have not found k2 NN b/c we have disconnected components at
        % the largest scale, simply move up/down to add additional
        % neighbors at infinite distance
        while (iup>1 && counter<k2+1)
            iup = iup-1;
            if counter==k2
                knndistances(s)=Inf;
            end
            counter=counter+1;
        end
        while (idown<nn && counter<k2+1)
            idown = idown+1;
            if counter==k2
                knndistances(s)=Inf;
            end
            counter=counter+1;
        end
    end

    N = size(sortedCCmatrix,1); UniqueDDs = unique(knndistances);
    SortedSimplices = sortedCCmatrix(:,m+1:end); 
    [~, idx] = unique(SortedSimplices, "rows"); Ths = 1.01*UniqueDDs; 
    Gknn = sparse(I, J, W, N, N); Gknn = max(Gknn, Gknn');

    tau = max(tau, 0.010);

    if d == 1 || d == 2 || d == 3 || d==4

        for i = 1:length(Ths) 
            Cutoff = Ths(i); 
            RowsToKeep = knndistances < Cutoff; % This corresponds to sorted simplices.
            RowsToKeep1 = RowsToKeep(idx);
            DenoisedGknn = Gknn(RowsToKeep1, RowsToKeep1); 
            [bins,binsizes] = conncomp(graph(DenoisedGknn));

            if max(binsizes) / length(bins) > 0.7 && max(binsizes) / N > 0.001 %> 0.70 
                DenoisedCCmatrix = sortedCCmatrix(RowsToKeep, :); break
                %if length(unique(DenoisedCCmatrix(:,31:34))) >= 0.5*5000, break; end
            end
        end
    else
        %for i=1:length(Ths)
        %    if sum(knndistances < Ths(i)) / length(knndistances) > 0.25
        %        Cutoff = Ths(i); break;
        %    end
        %end
        Epsilon=0.34;
        KeepPercent = 0; ll=find(Ths>(atan(tau/Epsilon*2))); Cutoff_Pos = ll(1);
        while KeepPercent < 0.05
            Cutoff = Ths(Cutoff_Pos); % 3 or 4
            %Cutoff = max(min(Ths(find(Ths>tau*1.05))), Ths(3));
            RowsToKeep = knndistances < Cutoff;
            KeepPercent = sum(RowsToKeep)/length(RowsToKeep);
            Cutoff_Pos = Cutoff_Pos+1;
        end
        %while KeepPercent < 0.05
        %    Cutoff = Ths(Cutoff_Pos); % 3 or 4
            %Cutoff = max(min(Ths(find(Ths>tau*1.05))), Ths(3));
        %    RowsToKeep = knndistances < Cutoff;
        %    KeepPercent = sum(RowsToKeep)/length(RowsToKeep);
        %    Cutoff_Pos = Cutoff_Pos+1;
        %end
        DenoisedCCmatrix = sortedCCmatrix(RowsToKeep, :);
        RowsToKeep1 = RowsToKeep(idx); DenoisedGknn = Gknn(RowsToKeep1, RowsToKeep1);
    end

%Cutoff = min(Ths(find(Ths>tau*1.05)));
%RowsToKeep = knndistances < Cutoff;
%DenoisedCCmatrix = sortedCCmatrix(RowsToKeep, :);

%RowsToKeep1 = RowsToKeep(idx); DenoisedGknn = Gknn(RowsToKeep1, RowsToKeep1);
DenoisedSimplices = unique(DenoisedCCmatrix(:,(m+1):end), 'rows');
[I,J,W] = find(DenoisedGknn); clear Gknn
[DenoisedCCmatrix, NewTh] = ConnectedComponents(I,J,W,DenoisedSimplices,NumScales);

end