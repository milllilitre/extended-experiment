% this file specifies the function used in least squares curve fit in
% processData.m
% x is real location and phase bias
% xData is phase vec after processing the measured data
% y is phase vec calculated

function yData = method2Func(x, xData)
tagPopulation = size(xData, 1) - 1;
tagLoc = xData(1:tagPopulation, :);
tags_each_array = tagPopulation / 2;

loc = [transpose(x(1:2)) xData((tagPopulation + 1), 3)];

channel = 920.625e6;
waveLength = 3e8 / channel;

Ddist = zeros(1, tagPopulation - 2);
Dphase = zeros(1, tagPopulation - 2);
for tagIndex = 1:(tags_each_array - 1)
    loc1 = tagLoc(tagIndex,:);
    loc2 = tagLoc(tagIndex + 1,:);
    dist1 = norm(loc1 - loc);
    dist2 = norm(loc2 - loc);
    Ddist(tagIndex) = dist2 - dist1;
    Dphase(tagIndex) = Ddist(tagIndex) / waveLength * 2 * pi;
end
for tagIndex = tags_each_array:(tagPopulation - 2)
    loc1 = tagLoc(tagIndex + 1,:);
    loc2 = tagLoc(tagIndex + 2,:);
    dist1 = norm(loc1 - loc);
    dist2 = norm(loc2 - loc);
    Ddist(tagIndex) = dist2 - dist1;
    Dphase(tagIndex) = Ddist(tagIndex) / waveLength * 2 * pi;
end
yData = transpose(Dphase);
end
