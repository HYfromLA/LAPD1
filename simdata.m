function [Data, LabelsGT] = simdata(DATAopts, j)

rng(j)  %rng('default')

%% 1d Examples

if strcmp(DATAopts.Shape, 'Two Lines')
    n1=DATAopts.Number(1); n2=DATAopts.Number(2); 
    D = DATAopts.AmbDim; Sigma = DATAopts.NoiseSigma; Theta = DATAopts.Angles;

    Rotation1 = eye(D, D); Rotation1(1,1) = cos(Theta(1));  Rotation1(1,2) = sin(Theta(1));
    Rotation1(2,1) = -sin(Theta(1)); Rotation1(2,2) = cos(Theta(1));
    Rotation2 = eye(D, D); Rotation2(1,1) = cos(Theta(2));  Rotation2(1,2) = sin(Theta(2));
    Rotation2(2,1) = -sin(Theta(2)); Rotation2(2,2) = cos(Theta(2));
    Effective1 = -.5+1*rand(n1,1); Noise1 = normrnd(0,Sigma, [n1, D]); %Noise1 = -Sigma + rand(n1,D)*2*Sigma; %Noise1 = normrnd(0,Sigma, [n1, D]);
    Data1 = Noise1; Data1(:, 1) = Data1(:, 1) + Effective1;
    Effective2 = -.5+1*rand(n2,1); Noise2 = normrnd(0,Sigma, [n2, D]); %Noise2 = -Sigma + rand(n2,D)*2*Sigma; %Noise2 = normrnd(0,Sigma, [n2, D]);
    Data2 = Noise2; Data2(:, 1) = Data2(:, 1) + Effective2;
    Data1 = Data1 * Rotation1; Data2 = Data2 * Rotation2; 
    Data = cat(1, Data1, Data2); LabelsGT = [repelem(1, n1) repelem(2, n2)]';
    %if D == 2 
    %    figure; plot(Data(1:n1,1),Data(1:n1, 2), '.', Data(n1+1:n1+n2,1),Data(n1+1:n1+n2,2),'.');
    %    axis equal; axis off; xlim([-.7, .7]); ylim([-.7, .7])
    %end
end

if strcmp(DATAopts.Shape, 'Two Curves')
    n1=DATAopts.Number(1); n2=DATAopts.Number(2); D = DATAopts.AmbDim; 
    Sigma = DATAopts.NoiseSigma; Theta = DATAopts.Angles; R = 1/DATAopts.Curvature;
    
    Rotation1 = eye(D, D); Rotation1(1,1) = cos(Theta(1));  Rotation1(1,2) = sin(Theta(1));
    Rotation1(2,1) = -sin(Theta(1)); Rotation1(2,2) = cos(Theta(1));
    Rotation2 = eye(D, D); Rotation2(1,1) = cos(Theta(2));  Rotation2(1,2) = sin(Theta(2));
    Rotation2(2,1) = -sin(Theta(2)); Rotation2(2,2) = cos(Theta(2));
    
    %tts = (pi/2)-asin(0.5/R) + 2*asin(0.5/R)*rand(n1, 2);
    tts = (pi/2)-0.5/R + 2*0.5/R*rand(n1, 2);
    %tt1=(pi/2)*rand(n1,1)+pi/4; tt2=(pi/2)*rand(n2,1)+pi/4;
    Effective1 = [R*cos(tts(:,1)) R*sin(tts(:,1))];  Noise1 = normrnd(0,Sigma, [n1, D]);
    Data1 = Noise1; Data1(:, [1 2]) = Data1(:, [1 2]) + Effective1;
    Effective2 = [R*cos(tts(:,2)) R*sin(tts(:,2))];  Noise2 = normrnd(0,Sigma, [n2, D]);
    Data2 = Noise2; Data2(:, [1 2]) = Data2(:, [1 2]) + Effective2;
    Data1 = Data1 * Rotation1; Data2 = Data2 * Rotation2; 
    Data2(:,1) = Data2(:,1)+R*sin(Theta(2)); Data2(:,2) = Data2(:,2)+R-R*cos(Theta(2)); 
    Data = cat(1, Data1, Data2); LabelsGT = [repelem(1, n1) repelem(2, n2)]';
    if D == 2
        figure; plot(Data(1:n1,1),Data(1:n1,2), '.', Data(n1+1:n1+n2,1),Data(n1+1:n1+n2,2),'.');
        axis equal; xlim([-.1 .2]); ylim([-.1 .2]);
        axis off; %xlim([-.1 .2]); ylim([-.1 .2]);
    end
end

if strcmp(DATAopts.Shape, 'Dollar Sign')
    % n1 is the S, n2 is the |. 
    n1 = DATAopts.Number(1); n2 = DATAopts.Number(2); D = DATAopts.AmbDim; Sigma = DATAopts.NoiseSigma;
    angle = 1.5*rand(n1/2, 2)*pi; X11=2*[cos(angle(:,1)) sin(angle(:,1))+1];
    X12=2*[-cos(angle(:,2)) -sin(angle(:,2))-1];
    X2 =2*[repelem(0,n2)' 5*rand(n2, 1)-2.5]; X = [X11; X12; X2];
    Noise = normrnd(0,Sigma, [n1+n2, D]); Noise(:, [1 2]) = Noise(:, [1 2])+X;
    Data=Noise; LabelsGT = [repelem(1, n1) repelem(2, n2)]';
    if D == 2
       figure; plot(Data(1:n1,1), Data(1:n1,2),'.', Data(n1+1:n1+n2,1), Data(n1+1:n1+n2,2),'.')
       axis equal; axis off
    end
end

%{
if strcmp(DATAopts.Shape, 'Dollar Sign')
     N=10000; D = DATAopts.AmbDim; Sigma = DATAopts.NoiseSigma;
     p = (3*pi)/(6+3*pi); X = zeros(N,2); Y = zeros(N,1); m = zeros(N,1); 
     for i = 1:N
        if rand() > p
            m(i,1) = 0; X(i,:) = [0 5*rand()-2.5];
            %Y(i,1) = 13 + randn(); % will be from 10 to 16 +/- 2
        else
            m(i,1) = 1; angle = 1.5*rand()*pi;
            x = cos(angle); y = sin(angle);
	
            % choose between top and bottom of S
            if rand() > 0.5
                X(i,:) = [x y+1]; %Y(i,1) = angle;        % Y is from 0 to 3*pi
            else            
                X(i,:) = [-x -y-1]; %Y(i,1) = 3*pi - angle;
            end        
        end
     end

     tt1 = X(m==1, :); tt2 = X(m==0, :); n1 = length(m(m==1)); n2 = length(m(m==0));
     PreData1 = normrnd(0,Sigma, [n1, D]); PreData2 = normrnd(0,Sigma, [n2, D]); 
     PreData1(:, [1 2]) = PreData1(:, [1 2])+tt1; PreData2(:, [1 2]) = PreData2(:, [1 2])+tt2; 
     Data=[PreData1;PreData2]; LabelsGT = [repelem(1, n1) repelem(2, n2)]';
     if D == 2
         figure; plot(PreData1(:,1), PreData1(:,2),'.', PreData2(:,1), PreData2(:,2),'.')
         axis equal; axis off
     end
end
%}

if strcmp(DATAopts.Shape, 'Olympic Rings')
    n1=DATAopts.Number(1); n2=DATAopts.Number(2); n3=DATAopts.Number(3); n4=DATAopts.Number(4); n5=DATAopts.Number(5);
    D = DATAopts.AmbDim; Sigma = DATAopts.NoiseSigma; Angle = linspace(pi/4,9*pi/4,n1);
    %thetab = -pi/4; thetak = -pi/4; thetar = -pi/4;

    xb = 2*(cos(Angle) * 1); yb = 2*(sin(Angle) * 1); %zb = cos(Angle + thetab) * 0.1;

    xy = 2*(cos(Angle) * 1 + 10/9); yy = 2*(sin(Angle) * 1 - 10/9);

    xk = 2*(cos(Angle) * 1 + 20/9); yk = 2*(sin(Angle) * 1);  %zk = cos(Angle + thetak) * 0.1;

    xg = 2*(cos(Angle) * 1 + 30/9); yg = 2*(sin(Angle) * 1 - 10/9);

    xr = 2*(cos(Angle) * 1 + 40/9); yr = 2*(sin(Angle) * 1); %zr = cos(Angle + thetar) * 0.1;

    PreData1 = normrnd(0,Sigma, [n1, D]); PreData2 = normrnd(0,Sigma, [n2, D]);
    PreData3 = normrnd(0,Sigma, [n3, D]); PreData4 = normrnd(0,Sigma, [n4, D]);
    PreData5 = normrnd(0,Sigma, [n5, D]); 
    PreData1(:,[1 2]) = PreData1(:,[1 2])+[xb',yb']; PreData2(:,[1 2]) = PreData2(:,[1 2])+[xy',yy']; 
    PreData3(:,[1 2]) = PreData3(:,[1 2])+[xk',yk']; PreData4(:,[1 2]) = PreData4(:,[1 2])+[xg',yg']; 
    PreData5(:,[1 2]) = PreData5(:,[1 2])+[xr',yr']; 
    Data = [PreData1; PreData2; PreData3; PreData4; PreData5]; LabelsGT = [repelem(1, n1) repelem(2, n2) repelem(3, n3) repelem(4, n4) repelem(5, n5)]';

    if D == 2
        h2 = figure;
        hold on
        plot(PreData1(:,1),PreData1(:,2),'b.','markersize',10);
        plot(PreData2(:,1),PreData2(:,2),'y.','markersize',10);
        plot(PreData3(:,1),PreData3(:,2),'k.','markersize',10);
        plot(PreData4(:,1),PreData4(:,2),'g.','markersize',10);
        plot(PreData5(:,1),PreData5(:,2),'r.','markersize',10);
        % make the axis pretty
        axis equal
        axis off
        %xlim([-1.2 5.2])
        set(h2,'Color',[1 1 1])
        hold off
    end
end

if strcmp(DATAopts.Shape, 'Three Curves')
    n1=DATAopts.Number(1); n2=DATAopts.Number(2); n3=DATAopts.Number(3);
    Sigma = DATAopts.NoiseSigma; Theta = DATAopts.Angles; D = DATAopts.AmbDim; 
    Rotation1 = eye(D, D); Rotation1(1,1) = cos(Theta(1));  Rotation1(1,2) = sin(Theta(1));
    Rotation1(2,1) = -sin(Theta(1)); Rotation1(2,2) = cos(Theta(1));
    Rotation2 = eye(D, D); Rotation2(1,1) = cos(Theta(2));  Rotation2(1,2) = sin(Theta(2));
    Rotation2(2,1) = -sin(Theta(2)); Rotation2(2,2) = cos(Theta(2));
    Rotation3 = eye(D, D); Rotation3(1,1) = cos(Theta(3));  Rotation3(1,2) = sin(Theta(3));
    Rotation3(2,1) = -sin(Theta(3)); Rotation3(2,2) = cos(Theta(3));

    r = 2; 

    theta1 = -pi+pi*rand(n1,1); x_1 = r*cos(theta1); z_1 = r*sin(theta1)+1; PreData1 = [x_1 z_1]; 
    theta2 = -pi+pi*rand(n2,1); x_2 = r*cos(theta2); z_2 = r*sin(theta2)+1; PreData2 = [x_2 z_2]; 
    theta3 = -pi+pi*rand(n3,1); x_3 = r*cos(theta3); z_3 = r*sin(theta3)+1; PreData3 = [x_3 z_3]; 
    Data1 = normrnd(0,Sigma, [n1, D]); Data2 = normrnd(0,Sigma, [n2, D]); Data3 = normrnd(0,Sigma, [n3, D]);
    Data1(:, [1 3]) = Data1(:, [1 3]) + PreData1; Data1 = Data1 * Rotation1;
    Data2(:, [1 3]) = Data2(:, [1 3]) + PreData2; Data2 = Data2 * Rotation2;
    Data3(:, [1 3]) = Data3(:, [1 3]) + PreData3; Data3 = Data3 * Rotation3;

    Data = cat(1, Data1, Data2, Data3); 
    LabelsGT = [repelem(1, n1) repelem(2, n2) repelem(3, n3)]';
    if D == 3
        figure; plot3(Data(1:n1,1),Data(1:n1,2),Data(1:n1,3),'.', ...
            Data(n1+1:n1+n2,1),Data(n1+1:n1+n2,2),Data(n1+1:n1+n2,3),'.', ...
            Data(n1+n2+1:n1+n2+n3,1),Data(n1+n2+1:n1+n2+n3,2),Data(n1+n2+1:n1+n2+n3,3),'.');
        axis equal; axis off
    end
end

%% 2d Examples

if strcmp(DATAopts.Shape, 'Two Curved Surfaces')
    n1=DATAopts.Number(1); n2=DATAopts.Number(2); D = DATAopts.AmbDim; 
    Sigma = DATAopts.NoiseSigma; Theta = DATAopts.Angles; R = 1/DATAopts.Curvature;
    
    Rotation1 = eye(D, D); Rotation1(1,1) = cos(Theta(1));  Rotation1(1,3) = sin(Theta(1));
    Rotation1(3,1) = -sin(Theta(1)); Rotation1(3,3) = cos(Theta(1));
    Rotation2 = eye(D, D); Rotation2(1,1) = cos(Theta(2));  Rotation2(1,3) = sin(Theta(2));
    Rotation2(3,1) = -sin(Theta(2)); Rotation2(3,3) = cos(Theta(2));
    
    %tts = (pi/2)-asin(0.5/R) + 2*asin(0.5/R)*rand(n1, 2);
    tts = (pi/2)-0.5/R + 2*0.5/R*rand(n1, 2);
    %tt1=(pi/2)*rand(n1,1)+pi/4; tt2=(pi/2)*rand(n2,1)+pi/4;
    Effective1 = [R*cos(tts(:,1)) -.25+.5*rand(n1,1) R*sin(tts(:,1))];  Noise1 = normrnd(0,Sigma, [n1, D]);
    Data1 = Noise1; Data1(:, [1 2 3]) = Data1(:, [1 2 3]) + Effective1;
    Effective2 = [R*cos(tts(:,2)) -.25+.5*rand(n2,1) R*sin(tts(:,2))];  Noise2 = normrnd(0,Sigma, [n2, D]);
    Data2 = Noise2; Data2(:, [1 2 3]) = Data2(:, [1 2 3]) + Effective2;
    Data1 = Data1 * Rotation1; Data2 = Data2 * Rotation2; 
    Data2(:,1) = Data2(:,1)+R*sin(Theta(2)); Data2(:,3) = Data2(:,3)+(R-R*cos(Theta(2))); 
    Data = cat(1, Data1, Data2); LabelsGT = [repelem(1, n1) repelem(2, n2)]';
    if D == 3
        figure; plot3(Data(1:n1,1),Data(1:n1,2),Data(1:n1,3), '.', Data(n1+1:n1+n2,1),Data(n1+1:n1+n2,2),Data(n1+1:n1+n2,3),'.');
        axis equal; axis off; xlim([-.15, .2]); zlim([-.15, .2]);
    end
end

if strcmp(DATAopts.Shape, 'Two Curved Cubes')
    n1=DATAopts.Number(1); n2=DATAopts.Number(2); D = DATAopts.AmbDim; 
    Sigma = DATAopts.NoiseSigma; Theta = DATAopts.Angles; R = 1/DATAopts.Curvature;
    
    Rotation1 = eye(D, D); Rotation1(1,1) = cos(Theta(1));  Rotation1(1,4) = sin(Theta(1));
    Rotation1(4,1) = -sin(Theta(1)); Rotation1(4,4) = cos(Theta(1));
    Rotation2 = eye(D, D); Rotation2(1,1) = cos(Theta(2));  Rotation2(1,4) = sin(Theta(2));
    Rotation2(4,1) = -sin(Theta(2)); Rotation2(4,4) = cos(Theta(2));
    
    %tts = (pi/2)-asin(0.5/R) + 2*asin(0.5/R)*rand(n1, 2);
    tts = (pi/2)-0.5/R + 2*0.5/R*rand(n1, 2);
    %tt1=(pi/2)*rand(n1,1)+pi/4; tt2=(pi/2)*rand(n2,1)+pi/4;
    Effective1 = [R*cos(tts(:,1)) -.25+.5*rand(n1,1) -.25+.5*rand(n1,1) R*sin(tts(:,1))];  Noise1 = normrnd(0,Sigma, [n1, D]);
    Data1 = Noise1; Data1(:, [1 2 3 4]) = Data1(:, [1 2 3 4]) + Effective1;
    Effective2 = [R*cos(tts(:,2)) -.25+.5*rand(n2,1) -.25+.5*rand(n1,1) R*sin(tts(:,2))];  Noise2 = normrnd(0,Sigma, [n2, D]);
    Data2 = Noise2; Data2(:, [1 2 3 4]) = Data2(:, [1 2 3 4]) + Effective2;
    Data1 = Data1 * Rotation1; Data2 = Data2 * Rotation2; 
    Data2(:,1) = Data2(:,1)+R*sin(Theta(2)); Data2(:,4) = Data2(:,4)+R-R*cos(Theta(2)); 
    Data = cat(1, Data1, Data2); LabelsGT = [repelem(1, n1) repelem(2, n2)]';
end

if strcmp(DATAopts.Shape, 'Two Planes')
    n1=DATAopts.Number(1); n2=DATAopts.Number(2); D = DATAopts.AmbDim; Sigma = DATAopts.NoiseSigma; Theta = DATAopts.Angles;
    Rotation1 = eye(D, D); Rotation1(1,1) = cos(Theta(1));  Rotation1(1,3) = sin(Theta(1));
    Rotation1(3,1) = -sin(Theta(1)); Rotation1(3,3) = cos(Theta(1));
    Rotation2 = eye(D, D); Rotation2(1,1) = cos(Theta(2));  Rotation2(1,3) = sin(Theta(2));
    Rotation2(3,1) = -sin(Theta(2)); Rotation2(3,3) = cos(Theta(2));
    Effective1 = -.5+1*rand(n1,2); Noise1 = normrnd(0,Sigma, [n1, D]); %-Sigma + rand(n1,D)*2*Sigma;
    Data1 = Noise1; Data1(:, [1 2]) = Data1(:, [1 2]) + Effective1;
    Effective2 = -.5+1*rand(n2,2); Noise2 = normrnd(0,Sigma, [n2, D]); %-Sigma + rand(n2,D)*2*Sigma;
    Data2 = Noise2; Data2(:, [1 2]) = Data2(:, [1 2]) + Effective2;
    Data1 = Data1 * Rotation1; Data2 = Data2 * Rotation2; 
    Data = cat(1, Data1, Data2); LabelsGT = [repelem(1, n1) repelem(2, n2)]';

    %if D == 3
    %    figure; plot3(Data(1:n1,1),Data(1:n1,2),Data(1:n1,3), '.', Data(n1+1:n1+n2,1),Data(n1+1:n1+n2,2),Data(n1+1:n1+n2,3),'.');
    %    xlabel('x'); ylabel('y'); zlabel('z'); axis equal; axis off; 
    %    xlim([-.7, .7]);ylim([-.7, .7]);zlim([-.7, .7]);
    %end
end


if strcmp(DATAopts.Shape, 'Two Planes Scaled')
    n1=DATAopts.Number(1); n2=DATAopts.Number(2); D = DATAopts.AmbDim; Sigma = DATAopts.NoiseSigma; Theta = DATAopts.Angles;
    Rotation1 = eye(D, D); Rotation1(1,1) = cos(Theta(1));  Rotation1(1,3) = sin(Theta(1));
    Rotation1(3,1) = -sin(Theta(1)); Rotation1(3,3) = cos(Theta(1));
    Rotation2 = eye(D, D); Rotation2(1,1) = cos(Theta(2));  Rotation2(1,3) = sin(Theta(2));
    Rotation2(3,1) = -sin(Theta(2)); Rotation2(3,3) = cos(Theta(2));
    Effective1 = -1+2*rand(n1,2); Noise1 = normrnd(0,Sigma, [n1, D]); %-Sigma + rand(n1,D)*2*Sigma;
    Data1 = Noise1; Data1(:, [1 2]) = Data1(:, [1 2]) + Effective1;
    Effective2 = -1+2*rand(n2,2); Noise2 = normrnd(0,Sigma, [n2, D]); %-Sigma + rand(n2,D)*2*Sigma;
    Data2 = Noise2; Data2(:, [1 2]) = Data2(:, [1 2]) + Effective2;
    Data1 = Data1 * Rotation1; Data2 = Data2 * Rotation2; 
    Data = cat(1, Data1, Data2); LabelsGT = [repelem(1, n1) repelem(2, n2)]';

    %if D == 3
    %    figure; plot3(Data(1:n1,1),Data(1:n1,2),Data(1:n1,3), '.', Data(n1+1:n1+n2,1),Data(n1+1:n1+n2,2),Data(n1+1:n1+n2,3),'.');
    %    xlabel('x'); ylabel('y'); zlabel('z'); axis equal; %axis off; 
    %end
end


%{
if strcmp(DATAopts.Shape, 'Two Triangles') % Testing intersection is of dimension 0. 
    n1=DATAopts.Number(1); n2=DATAopts.Number(2); D = DATAopts.AmbDim; Sigma = DATAopts.NoiseSigma; Theta = DATAopts.Angles;
    Rotation1 = eye(D, D); Rotation1(1,1) = cos(Theta(1));  Rotation1(1,3) = sin(Theta(1));
    Rotation1(3,1) = -sin(Theta(1)); Rotation1(3,3) = cos(Theta(1));
    Rotation2 = eye(D, D); Rotation2(1,1) = cos(Theta(2));  Rotation2(1,3) = sin(Theta(2));
    Rotation2(3,1) = -sin(Theta(2)); Rotation2(3,3) = cos(Theta(2));
    
    Effective1 = 0.5*rand(n1,1); Effective1(:,2) = 0; 
    for i=1:n1
        temp = Effective1(i, 1)*tan(pi/6)+Sigma;
        Effective1(i,2) = -temp + 2*temp * rand(1,1);
    end    
    Noise1 = normrnd(0,Sigma, [n1, D]);
    Data1 = Noise1; Data1(:, [1 2]) = Data1(:, [1 2]) + Effective1;

    Effective2 = 0.5*rand(n2,1); Effective2(:,2) = 0; 
    for i=1:n2
        temp = Effective2(i, 1)*tan(pi/6);
        Effective2(i,2) = -temp + 2*temp * rand(1,1);
    end    
    Noise2 = normrnd(0,Sigma, [n2, D]);
    Data2 = Noise2; Data2(:, [1 2]) = Data2(:, [1 2]) + Effective2;
    Data1 = Data1 * Rotation1; Data2 = Data2 * Rotation2; 
    Data = cat(1, Data1, Data2); LabelsGT = [repelem(1, n1) repelem(2, n2)]';

    if D == 3
        figure; plot3(Data(1:n1,1),Data(1:n1,2),Data(1:n1,3), '.', Data(n1+1:n1+n2,1),Data(n1+1:n1+n2,2),Data(n1+1:n1+n2,3),'.');
        xlabel('x')
        ylabel('y')  
        zlabel('z') 
        axis equal; %axis off; 
    end
end
%}


if strcmp(DATAopts.Shape, 'Two Triangles') % Testing intersection is of dimension 0. 
    n1=DATAopts.Number(1); n2=DATAopts.Number(2); D = DATAopts.AmbDim; Sigma = DATAopts.NoiseSigma; Theta = DATAopts.Angles;
    Rotation1 = eye(D, D); Rotation1(1,1) = cos(Theta(1));  Rotation1(1,3) = sin(Theta(1));
    Rotation1(3,1) = -sin(Theta(1)); Rotation1(3,3) = cos(Theta(1));
    Rotation2 = eye(D, D); Rotation2(1,1) = cos(Theta(2));  Rotation2(1,3) = sin(Theta(2));
    Rotation2(3,1) = -sin(Theta(2)); Rotation2(3,3) = cos(Theta(2));
    
    Effective1(:,1) = 1.6*(0.5*rand(n1,1)); Effective1(:,2) =1.6*(-.5+1*rand(n1,1));
    Effective1 = Effective1(abs(Effective1(:,2)) <= Effective1(:,1),:);
    n1 = size(Effective1, 1);
    Noise1 = normrnd(0,Sigma, [n1, D]);
    Data1 = Noise1; Data1(:, [1 2]) = Data1(:, [1 2]) + Effective1;

    Effective2(:,1) = 1.6*(0.5*rand(n2,1)); Effective2(:,2) = 1.6*(-.5+1*rand(n2,1));
    Effective2 = Effective2(abs(Effective2(:,2)) <= Effective2(:,1),:);
    n2 = size(Effective2, 1); 
    Noise2 = normrnd(0,Sigma, [n2, D]);
    Data2 = Noise2; Data2(:, [1 2]) = Data2(:, [1 2]) + Effective2;
    
    Data1 = Data1 * Rotation1; Data2 = Data2 * Rotation2; 
    Data = cat(1, Data1, Data2); LabelsGT = [repelem(1, n1) repelem(2, n2)]';

    if D == 3
        figure; plot3(Data(1:n1,1),Data(1:n1,2),Data(1:n1,3), '.', Data(n1+1:n1+n2,1),Data(n1+1:n1+n2,2),Data(n1+1:n1+n2,3),'.');
        %xlabel('x'), ylabel('y'), zlabel('z') 
        axis equal; axis off; 
    end
end



if strcmp(DATAopts.Shape, 'Two Spheres')
    n1 = DATAopts.Number(1); n2 = DATAopts.Number(2); D = DATAopts.AmbDim; Sigma = DATAopts.NoiseSigma;
    
    u=rand(n1,1); v=rand(n1,1); r=1; phi=2*pi*u; teta=acos(2*v-1);
    z=r*cos(teta); x=sqrt(r^2-z.^2).*cos(phi); y=sqrt(r^2-z.^2).*sin(phi);
    D1=[x y z];

    u=rand(n2,1); v=rand(n2,1); r=1; phi=2*pi*u; teta=acos(2*v-1);
    z=r*cos(teta); x=sqrt(r^2-z.^2).*cos(phi); y=sqrt(r^2-z.^2).*sin(phi)+1.5;
    D2=[x y z];

    PreData1 = normrnd(0,Sigma, [n1, D]); PreData2 = normrnd(0,Sigma, [n1, D]); 
    PreData1(:, [1 2 3]) = PreData1(:, [1 2 3])+D1; PreData2(:, [1 2 3]) = PreData2(:, [1 2 3])+D2; 
    Data=[PreData1;PreData2]; LabelsGT = [repelem(1, n1) repelem(2, n1)]';
    if D == 3
        figure; plot3(Data(1:n1,1), Data(1:n1,2), Data(1:n1,3), '.', Data(n1+1:n1+n2,1), Data(n1+1:n1+n2,2), Data(n1+1:n1+n2,3), '.')
        axis equal; axis off
    end
end

if strcmp(DATAopts.Shape, 'Swiss Roll')
    n1=DATAopts.Number(1); n2=DATAopts.Number(2); D = DATAopts.AmbDim; Sigma = DATAopts.NoiseSigma;
    D1(:,1)=-2.5+5*rand(n1,1); r =4*rand(n1,1); D1(:,2)=cos(0.8*pi*r).*r; D1(:,3)=-sin(0.8*pi*r).*r; 
    D2(:,1)=-2.5+5*rand(n2,1); D2(:,2)=-5+10.*rand(n2,1); 
    EffCols1 = D1; EffCols2 = D2; PreData1 = normrnd(0,Sigma, [n1, D]); PreData2 = normrnd(0,Sigma, [n2, D]); 
    PreData1(:,[1 2 3]) = PreData1(:,[1 2 3])+EffCols1; PreData2(:,[1 2]) = PreData2(:,[1 2])+EffCols2;
   
    Data = cat(1, PreData1, PreData2); LabelsGT = [repelem(1, n1) repelem(2, n2)]';

    if D == 3
        figure; plot3(Data(1:n1,1),Data(1:n1,2),Data(1:n1,3), '.', Data(n1+1:n1+n2,1),Data(n1+1:n1+n2,2),Data(n1+1:n1+n2,3),'.');
        %xlabel('x'), ylabel('y'), zlabel('z'), title('2-d Noisy Swiss Roll Demo', 'FontSize',16) 
        axis equal; axis off; 
    end
end

if strcmp(DATAopts.Shape, 'Three Planes(1)') % Three planes intersect at the same point. 
    n1=DATAopts.Number(1); n2=DATAopts.Number(2); n3=DATAopts.Number(3);
    D = DATAopts.AmbDim; Sigma = DATAopts.NoiseSigma; Theta = DATAopts.Angles;

    Rotation1 = eye(D, D); Rotation1(1,1) = cos(Theta(1));  Rotation1(1,3) = sin(Theta(1));
    Rotation1(3,1) = -sin(Theta(1)); Rotation1(3,3) = cos(Theta(1));
    Rotation2 = eye(D, D); Rotation2(1,1) = cos(Theta(2));  Rotation2(1,3) = sin(Theta(2));
    Rotation2(3,1) = -sin(Theta(2)); Rotation2(3,3) = cos(Theta(2));
    Rotation3 = eye(D, D); Rotation3(1,1) = cos(Theta(3));  Rotation3(1,3) = sin(Theta(3));
    Rotation3(3,1) = -sin(Theta(3)); Rotation3(3,3) = cos(Theta(3));

    Effective1 = -.5+1*rand(n1,2); Noise1 = normrnd(0,Sigma, [n1, D]);
    Data1 = Noise1; Data1(:, [1 2]) = Data1(:, [1 2]) + Effective1;
    Effective2 = -.5+1*rand(n2,2); Noise2 = normrnd(0,Sigma, [n2, D]);
    Data2 = Noise2; Data2(:, [1 2]) = Data2(:, [1 2]) + Effective2;
    Effective3 = -.5+1*rand(n3,2); Noise3 = normrnd(0,Sigma, [n3, D]);
    Data3 = Noise3; Data3(:, [1 2]) = Data3(:, [1 2]) + Effective3;
    Data1 = Data1 * Rotation1; Data2 = Data2 * Rotation2; Data3 = Data3 * Rotation3;
    Data = cat(1, Data1, Data2, Data3); LabelsGT = [repelem(1, n1) repelem(2, n2) repelem(3, n3)]';

    if D == 3
        figure; plot3(Data(1:n1,1),Data(1:n1,2),Data(1:n1,3), '.', Data(n1+1:n1+n2,1),Data(n1+1:n1+n2,2),Data(n1+1:n1+n2,3),'.', ...
            Data(n1+n2+1:n1+n2+n3,1),Data(n1+n2+1:n1+n2+n3,2),Data(n1+n2+1:n1+n2+n3,3),'.');
        %xlabel('x'), ylabel('y'), zlabel('z')
        axis equal, axis off
    end
end

if strcmp(DATAopts.Shape, 'Three Planes(2)') % Three planes intersect mutually. 
    n1=DATAopts.Number(1); n2=DATAopts.Number(2); n3=DATAopts.Number(3);
    D = DATAopts.AmbDim; Sigma = DATAopts.NoiseSigma; Theta = DATAopts.Angles;

    Rotation1 = eye(D, D); Rotation1(1,1) = cos(Theta(1));  Rotation1(1,3) = sin(Theta(1));
    Rotation1(3,1) = -sin(Theta(1)); Rotation1(3,3) = cos(Theta(1));
    Rotation2 = eye(D, D); Rotation2(1,1) = cos(Theta(2));  Rotation2(1,3) = sin(Theta(2));
    Rotation2(3,1) = -sin(Theta(2)); Rotation2(3,3) = cos(Theta(2));
    Rotation3 = eye(D, D); Rotation3(1,1) = cos(Theta(3));  Rotation3(1,3) = sin(Theta(3));
    Rotation3(3,1) = -sin(Theta(3)); Rotation3(3,3) = cos(Theta(3));

    Effective1 = -.5+1*rand(n1,2); Noise1 = normrnd(0,Sigma, [n1, D]);
    Data1 = Noise1; Data1(:, [1 2]) = Data1(:, [1 2]) + Effective1;
    Effective2 = -.5+1*rand(n2,2); Noise2 = normrnd(0,Sigma, [n2, D]);
    Data2 = Noise2; Data2(:, [1 2]) = Data2(:, [1 2]) + Effective2;
    Effective3 = -.5+1*rand(n3,2); Noise3 = normrnd(0,Sigma, [n3, D]);
    Data3 = Noise3; Data3(:, [1 2]) = Data3(:, [1 2]) + Effective3;
    Data1 = Data1 * Rotation1; Data2 = Data2 * Rotation2; Data3 = Data3 * Rotation3;
    Data2(:,1) = Data2(:,1) - 0.2; Data3(:,1) = Data3(:,1) + 0.2;
    Data2(:,3) = Data2(:,3) + 0.25; Data3(:,3) = Data3(:,3) + 0.25;
    Data = cat(1, Data1, Data2, Data3); LabelsGT = [repelem(1, n1) repelem(2, n2) repelem(3, n3)]';

    if D == 3
        figure; plot3(Data(1:n1,1),Data(1:n1,2),Data(1:n1,3), '.', Data(n1+1:n1+n2,1),Data(n1+1:n1+n2,2),Data(n1+1:n1+n2,3),'.', ...
            Data(n1+n2+1:n1+n2+n3,1),Data(n1+n2+1:n1+n2+n3,2),Data(n1+n2+1:n1+n2+n3,3),'.');
        axis equal; axis off
    end
end

if strcmp(DATAopts.Shape, 'Three Planes(3)') % Three planes intersect mutually. 
    n1=DATAopts.Number(1); n2=DATAopts.Number(2); n3=DATAopts.Number(3);
    D = DATAopts.AmbDim; Sigma = DATAopts.NoiseSigma; Theta = DATAopts.Angles;

    Rotation1 = eye(D, D); Rotation1(1,1) = cos(Theta(1));  Rotation1(1,3) = sin(Theta(1));
    Rotation1(3,1) = -sin(Theta(1)); Rotation1(3,3) = cos(Theta(1));
    Rotation2 = eye(D, D); Rotation2(1,1) = cos(Theta(2));  Rotation2(1,3) = sin(Theta(2));
    Rotation2(3,1) = -sin(Theta(2)); Rotation2(3,3) = cos(Theta(2));
    Rotation3 = eye(D, D); Rotation3(1,1) = cos(Theta(3));  Rotation3(1,3) = sin(Theta(3));
    Rotation3(3,1) = -sin(Theta(3)); Rotation3(3,3) = cos(Theta(3));

    Effective1(:,1) = 1*(-.75+1.5*rand(n1,1)); Effective1(:,2) = 1*(-.5+1*rand(n1,1)); Noise1 = normrnd(0,Sigma, [n1, D]);
    Data1 = Noise1; Data1(:, [1 2]) = Data1(:, [1 2]) + Effective1;
    Effective2 = 1*(-.5+1*rand(n2,2)); Noise2 = normrnd(0,Sigma, [n2, D]);
    Data2 = Noise2; Data2(:, [1 2]) = Data2(:, [1 2]) + Effective2;
    Effective3 = 1*(-.5+1*rand(n3,2)); Noise3 = normrnd(0,Sigma, [n3, D]);
    Data3 = Noise3; Data3(:, [1 2]) = Data3(:, [1 2]) + Effective3;
    Data1 = Data1 * Rotation1; Data2 = Data2 * Rotation2; Data3 = Data3 * Rotation3;
    Data2(:,1) = 1*(Data2(:,1) - 0.15); Data3(:,1) = 1*(Data3(:,1) + 0.15);
    Data2(:,3) = 1*(Data2(:,3) + 0.15); Data3(:,3) = 1*(Data3(:,3) + 0.15);
    Data = cat(1, Data1, Data2, Data3); LabelsGT = [repelem(1, n1) repelem(2, n2) repelem(3, n3)]';

    if D == 3
        figure; plot3(Data(1:n1,1),Data(1:n1,2),Data(1:n1,3), '.', Data(n1+1:n1+n2,1),Data(n1+1:n1+n2,2),Data(n1+1:n1+n2,3),'.', ...
            Data(n1+n2+1:n1+n2+n3,1),Data(n1+n2+1:n1+n2+n3,2),Data(n1+n2+1:n1+n2+n3,3),'.');
        axis equal; axis off
    end
end

%% 3d examples 
if strcmp(DATAopts.Shape, 'Two Cuboids')
    n1=DATAopts.Number(1); n2=DATAopts.Number(2); 
    D = DATAopts.AmbDim; Sigma = DATAopts.NoiseSigma; Theta = DATAopts.Angles;
    Rotation1 = eye(D, D); Rotation1(1,1) = cos(Theta(1));  Rotation1(1,4) = sin(Theta(1));
    Rotation1(4,1) = -sin(Theta(1)); Rotation1(4,4) = cos(Theta(1));
    Rotation2 = eye(D, D); Rotation2(1,1) = cos(Theta(2));  Rotation2(1,4) = sin(Theta(2));
    Rotation2(4,1) = -sin(Theta(2)); Rotation2(4,4) = cos(Theta(2));
    
    Effective1 = -.5+1*rand(n1,3); Noise1 = normrnd(0,Sigma, [n1, D]); %-Sigma + rand(n1,D)*2*Sigma;
    Data1 = Noise1; Data1(:, [1 2 3]) = Data1(:, [1 2 3]) + Effective1;
    Effective2 = -.5+1*rand(n2,3); Noise2 = normrnd(0,Sigma, [n2, D]); %-Sigma + rand(n2,D)*2*Sigma;
    Data2 = Noise2; Data2(:, [1 2 3]) = Data2(:, [1 2 3]) + Effective2;
    Data1 = Data1 * Rotation1; Data2 = Data2 * Rotation2; 
    Data = cat(1, Data1, Data2); LabelsGT = [repelem(1, n1) repelem(2, n2)]';
end

if strcmp(DATAopts.Shape, 'd=4')
    n1=DATAopts.Number(1); n2=DATAopts.Number(2); 
    D = DATAopts.AmbDim; Sigma = DATAopts.NoiseSigma; Theta = DATAopts.Angles;
    Rotation1 = eye(D, D); Rotation1(1,1) = cos(Theta(1));  Rotation1(1,5) = sin(Theta(1));
    Rotation1(5,1) = -sin(Theta(1)); Rotation1(5,5) = cos(Theta(1));
    Rotation2 = eye(D, D); Rotation2(1,1) = cos(Theta(2));  Rotation2(1,5) = sin(Theta(2));
    Rotation2(5,1) = -sin(Theta(2)); Rotation2(5,5) = cos(Theta(2));
    
    Effective1 = -.5+1*rand(n1,4); Noise1 = normrnd(0,Sigma, [n1, D]); %-Sigma + rand(n1,D)*2*Sigma;
    Data1 = Noise1; Data1(:, [1 2 3 4]) = Data1(:, [1 2 3 4]) + Effective1;
    Effective2 = -.5+1*rand(n2,4); Noise2 = normrnd(0,Sigma, [n2, D]); %-Sigma + rand(n2,D)*2*Sigma;
    Data2 = Noise2; Data2(:, [1 2 3 4]) = Data2(:, [1 2 3 4]) + Effective2;
    Data1 = Data1 * Rotation1; Data2 = Data2 * Rotation2; 
    Data = cat(1, Data1, Data2); LabelsGT = [repelem(1, n1) repelem(2, n2)]';
end

if strcmp(DATAopts.Shape, 'd=5')
    n1=DATAopts.Number(1); n2=DATAopts.Number(2); 
    D = DATAopts.AmbDim; Sigma = DATAopts.NoiseSigma; Theta = DATAopts.Angles;
    Rotation1 = eye(D, D); Rotation1(1,1) = cos(Theta(1));  Rotation1(1,6) = sin(Theta(1));
    Rotation1(6,1) = -sin(Theta(1)); Rotation1(6,6) = cos(Theta(1));
    Rotation2 = eye(D, D); Rotation2(1,1) = cos(Theta(2));  Rotation2(1,6) = sin(Theta(2));
    Rotation2(6,1) = -sin(Theta(2)); Rotation2(6,6) = cos(Theta(2));
    
    Effective1 = -.5+1*rand(n1,5); Noise1 = normrnd(0,Sigma, [n1, D]); %-Sigma + rand(n1,D)*2*Sigma;
    Data1 = Noise1; Data1(:, 1:5) = Data1(:, 1:5) + Effective1;
    Effective2 = -.5+1*rand(n2,5); Noise2 = normrnd(0,Sigma, [n2, D]); %-Sigma + rand(n2,D)*2*Sigma;
    Data2 = Noise2; Data2(:, 1:5) = Data2(:, 1:5) + Effective2;
    Data1 = Data1 * Rotation1; Data2 = Data2 * Rotation2; 
    Data = cat(1, Data1, Data2); LabelsGT = [repelem(1, n1) repelem(2, n2)]';
end


if strcmp(DATAopts.Shape, 'Rose & Circle')
    n1=DATAopts.Number(1); n2=DATAopts.Number(2); n3=DATAopts.Number(3);
    D=DATAopts.AmbDim; Sigma = DATAopts.NoiseSigma;
    teta=(2*pi)*randn(n1,1);     % Using  spherical coordinates
    x1=2*cos(teta/.5).*sin(teta);
    y1=2*cos(teta/.5).*cos(teta);
    D1=[x1 y1];

    teta=(2*pi)*randn(n2,1);     % Using  spherical coordinates
    x2=2*(.5*sin(teta));
    y2=2*(.5*cos(teta));
    D2=[x2 y2];

    teta=(2*pi)*randn(n3,1);     % Using  spherical coordinates
    x3=2*(.5*sin(teta)+.5);
    y3=2*(.5*cos(teta)+.5);
    D3=[x3 y3];

    PreData1 = normrnd(0,Sigma, [n1, D]); PreData2 = normrnd(0,Sigma, [n2, D]); PreData3 = normrnd(0,Sigma, [n3, D]);
    PreData1(:, [1 2]) = PreData1(:, [1 2]) + D1; PreData2(:, [1 2]) = PreData2(:,[1 2]) + D2; PreData3(:, [1 2]) = PreData3(:, [1 2])+D3; 
    Data=[PreData1;PreData2; PreData3]; LabelsGT = [repelem(1, n1) repelem(2, n2) repelem(3, n3)]';
    if D == 2
        figure; plot(PreData1(:,1),PreData1(:,2), '.', ...
            PreData2(:,1),PreData2(:,2),'.', PreData3(:,1),PreData3(:,2),'.');
        axis equal; axis off
    end
end

if strcmp(DATAopts.Shape, 'Three Lines')
    n1=DATAopts.Number(1); n2=DATAopts.Number(2); n3=DATAopts.Number(3);
    D = DATAopts.AmbDim; Sigma = DATAopts.NoiseSigma; Theta = DATAopts.Angles;

    Rotation1 = eye(D, D); Rotation1(1,1) = cos(Theta(1));  Rotation1(1,2) = sin(Theta(1));
    Rotation1(2,1) = -sin(Theta(1)); Rotation1(2,2) = cos(Theta(1));
    Rotation2 = eye(D, D); Rotation2(1,1) = cos(Theta(2));  Rotation2(1,2) = sin(Theta(2));
    Rotation2(2,1) = -sin(Theta(2)); Rotation2(2,2) = cos(Theta(2));
    Rotation3 = eye(D, D); Rotation3(1,1) = cos(Theta(3));  Rotation3(1,2) = sin(Theta(3));
    Rotation3(2,1) = -sin(Theta(3)); Rotation3(2,2) = cos(Theta(3));

    Effective1 = -.5+1*rand(n1,1); Noise1 = normrnd(0,Sigma, [n1, D]); %Noise1 = -Sigma + rand(n1,D)*2*Sigma; %Noise1 = normrnd(0,Sigma, [n1, D]);
    Data1 = Noise1; Data1(:, 1) = Data1(:, 1) + Effective1;
    Effective2 = -.5+1*rand(n2,1); Noise2 = normrnd(0,Sigma, [n2, D]); %Noise2 = -Sigma + rand(n2,D)*2*Sigma; %Noise2 = normrnd(0,Sigma, [n2, D]);
    Data2 = Noise2; Data2(:, 1) = Data2(:, 1) + Effective2;
    Effective3 = -.5+1*rand(n3,1); Noise3 = normrnd(0,Sigma, [n3, D]); %Noise2 = -Sigma + rand(n2,D)*2*Sigma; %Noise2 = normrnd(0,Sigma, [n2, D]);
    Data3 = Noise3; Data3(:, 1) = Data3(:, 1) + Effective3;
    Data1 = Data1 * Rotation1; Data2 = Data2 * Rotation2; Data3 = Data3 * Rotation3;
    Data = cat(1, Data1, Data2, Data3); LabelsGT = [repelem(1, n1) repelem(2, n2) repelem(3, n3)]';
    if D == 2 
        figure; plot(Data(1:n1,1),Data(1:n1, 2), '.', Data(n1+1:n1+n2,1),Data(n1+1:n1+n2,2),'.', ...
            Data(n1+n2+1:n1+n2+n3,1),Data(n1+n2+1:n1+n2+n3,2),'.');
        axis equal; axis off;
    end
end


end