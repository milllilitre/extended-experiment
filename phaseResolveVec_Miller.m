% phase based localization function library
% input and output are row vectors
% because miller based coding may cause reader's phase measurements to have
% a phase error of pi
function output = phaseResolveVec_Miller(input)
    % resolve phase, input is a row vector
    length = size(input,2);
    if(length > 0)
        output = input;
        delta = max(output) - min(output);
        if((delta > 0.9* pi) && (delta < 1.2 * pi))
            threshold = (max(output) + min(output)) / 2;
            tmpCount = sum(output > threshold);
            if tmpCount > length / 2
                for i = 1:length
                    if(output(i) < threshold)
                        output(i) = output(i) + pi;
                    end
                end
            else
                for i = 1:length
                    if(output(i) > threshold)
                        output(i) = output(i) - pi;
                    end
                end
            end
        end
    else
        output = 0;
    end
end
