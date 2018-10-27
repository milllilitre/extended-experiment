function [I,Q] = rssiPhase2IQ(rssi, phase)
amp = 10 .^ (rssi / 10);
I = amp .* cos(phase);
Q = amp .* sin(phase);