function [rssi, phase] = IQ2rssiPhase(I, Q)
amp = sqrt(I.^2 + Q.^2);
rssi = 10 .* log10(amp);
phase = mod(atan2(Q, I), 2 * pi);
end

