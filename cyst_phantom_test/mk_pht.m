%  Make the scatteres for a simulation and store
%  it in a file for later simulation use

%   Joergen Arendt Jensen, Feb. 26, 1997

[phantom_positions, phantom_amplitudes] = cyst_pht(1000);
save pht_data.mat phantom_positions phantom_amplitudes
