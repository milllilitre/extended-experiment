% phase based localization function library
% input and output are row vectors
function output = phaseResolveRefVec_Miller(input, reference)
%     figure();subplot(3,1,1);
%     plot(input);
    output = phaseResolveVec_Miller(input);
%     subplot(3,1,2);plot(output);
    tmpMean = mean(output);
    if(abs(tmpMean - pi - reference) < abs(tmpMean - reference))
        output = output - pi;
    end
    if(abs(tmpMean + pi - reference) < abs(tmpMean - reference))
        output = output + pi;
    end
end
