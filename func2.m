% this file specifies the function used in least squares curve fit in
% processData.m
% x is real location and phase bias
% xData is phase vec after processing the measured data
% y is phase vec calculated
function yData = func2(x, xData)
loc = [transpose(x(1:2)) 0.59];
tagPhase = xData;

tagPopulation = 34;
tagLoc = zeros(tagPopulation,3);
tagLoc(1:17,1) = 0.6:0.15:3;
tagLoc(1:17,2) = 3.6 * ones(17,1);
tagLoc(1:17,3) = ones(17,1) * 0.43;
tagLoc(18:34,1) = 0 * ones(17,1);
tagLoc(18:34,2) = 0.6:0.15:3;
tagLoc(18:34,3) = 0.43 * ones(17,1);

readerLoc = [1.8 0 0.59
    3.6 1.8 0.59];

channel = 920.625e6;
waveLength = 3e8 / channel;

Ddist = zeros(1, 16);
Dphase = zeros(1, 16);
phaseVec = zeros(1, 17);
phaseVec(1) = x(3);
for tagIndex = 1:16
    dist1 = norm(tagLoc(tagIndex + 17,:) - loc);
    dist2 = norm(tagLoc(tagIndex + 18,:) - loc);
    Ddist(tagIndex) = dist2 - dist1;
    Dphase(tagIndex) = Ddist(tagIndex) / waveLength * 2 * pi;
    phaseVec(tagIndex + 1) = phaseVec(tagIndex) + Dphase(tagIndex);
end
yData = transpose(phaseVec);
end
