% check upper and lower bounds of input, if a value is out of bound, then 
% set it to boundry
% input is a matrix of row vectors, upperBound and lowerBounds are row
% vectors, output is bounded matrix of row vectors
function [output] = checkBounds(input, upperBound, lowerBound)
output = input;
dataNum = size(output, 1);
upperBoundMat = repmat(upperBound, dataNum, 1);
lowerBoundMat = repmat(lowerBound, dataNum, 1);
tmpMat = output > upperBoundMat;
output(tmpMat) = upperBoundMat(tmpMat);
tmpMat = output < lowerBoundMat;
output(tmpMat) = lowerBoundMat(tmpMat);
end