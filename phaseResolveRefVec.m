% phase based localization function library
% input and output are row vectors
function output = phaseResolveRefVec(input, reference)
%     figure();subplot(3,1,1);
%     plot(input);
output = phaseResolveVec(input);
%     subplot(3,1,2);plot(output);
tmpMean = mean(output);
if(abs(tmpMean - 2 * pi - reference) < abs(tmpMean - reference))
    output = output - 2 * pi;
end
if(abs(tmpMean + 2 * pi - reference) < abs(tmpMean - reference))
    output = output + 2 * pi;
end
%     subplot(3,1,3);plot(output);
    
end
