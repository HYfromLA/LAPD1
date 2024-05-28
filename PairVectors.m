function [V1, V2] = PairVectors(d, m, temp)

if d == 1
    NormedVecs = normc(temp)';
    V1 = repelem(NormedVecs(1:(m-1), :), (m-1):-1:1, 1);
    V2 = repmat(NormedVecs, (m-1), 1);
else
    SharedVectors = temp(:,1:(d-1)); OtherVectors = temp(:, d:end);   % The shared vectors and the other vectors. 
    % Make shared vectors an orthonormal basis; use residuals from
    % projecting other vectors to the basis to calculate the angles. 
    P = orth(SharedVectors); Q = P';  
    Residuals = OtherVectors - P*(Q*OtherVectors); Residuals = Residuals ./ vecnorm(Residuals, 2, 1);    %Residuals = normc(OtherVectors - P*(Q*OtherVectors));
    Residuals = Residuals';
    V1 = repelem(Residuals(1:(m-1), :), (m-1):-1:1, 1);
    V2 = repmat(Residuals, (m-1), 1);
end

end