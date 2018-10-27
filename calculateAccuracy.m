function [distMat, avaDist, medDist, maxDist, variance] = calculateAccuracy(predictLoc, trueLoc)
diffMat = abs(trueLoc - predictLoc);
distMat = sqrt(sum(diffMat.^2, 2));
avaDist = nanmean(distMat);
medDist = nanmedian(distMat);
maxDist = nanmax(distMat);
tmpX = zeros(size(predictLoc));
tmpY = zeros(size(predictLoc));
tmpX(:,1) = predictLoc(:,1);
tmpX(:,2) = trueLoc(:,1);
tmpY(:,1) = predictLoc(:,2);
tmpY(:,2) = trueLoc(:,2);
variance = var(distMat);
% figure();hold on;
% tmpNum = size(predictLoc, 1);
%     for i = 1:tmpNum
%         plot(tmpX(i,:), tmpY(i,:));
%     end


end