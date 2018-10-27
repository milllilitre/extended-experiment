function [outRssi, outPhase] = vecSum(inRssi1, inPhase1, inRssi2, inPhase2)
    [inI1, inQ1] = rssiPhase2IQ(inRssi1, inPhase1);
    [inI2, inQ2] = rssiPhase2IQ(inRssi2, inPhase2);
    [outRssi, outPhase] = IQ2rssiPhase(inI1 + inI2, inQ1 + inQ2);
end




