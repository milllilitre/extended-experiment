% this file specifies the function used in least squares curve fit in
% processData.m
% x is real location and phase bias
% xData is phase vec after processing the measured data
% y is phase vec calculated
function yData = method2Func2(x, xData)
tagPopulation = size(xData, 1) - 1;
tagLoc = xData(1:tagPopulation, :);
tags_each_array = tagPopulation / 2;

loc = [transpose(x(1:2)) xData((tagPopulation + 1), 3)];

channel = 920.625e6;
waveLength = 3e8 / channel;

Ddist = zeros(1, (tags_each_array - 1));
Dphase = zeros(1, (tags_each_array - 1));
phaseVec = zeros(1, tags_each_array);
for tagIndex = 1:(tags_each_array - 1)
    dist1 = norm(tagLoc(tagIndex + tags_each_array,:) - loc);
    dist2 = norm(tagLoc(tagIndex + tags_each_array + 1,:) - loc);
    Ddist(tagIndex) = dist2 - dist1;
    Dphase(tagIndex) = Ddist(tagIndex) / waveLength * 2 * pi;
end
yData = transpose(Dphase);
end
