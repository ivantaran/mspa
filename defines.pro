
// mhz = 1.0e6;
speed_of_ligth = 299792458.0;
freq = GetNumber["freq"];
// freq = GetNumber["Model/Frequency"];
// freq = 140.0;
// freq = freq * mhz;
lambda = speed_of_ligth / freq;
k0 = 2.0 * Pi / lambda;
// epr = GetNumber["Model/epr"];
epr = 1.5;

SetNumber["Model/Lambda", lambda];
SetNumber["Model/WaveNumber", k0];
// Printf["ololo: %g", GetNumber["Model/Frequency"]];
Printf["ololo!"];
Printf[OnelabAction];

SKINFEED = 1;
SKINCONDUCTOR = 2;
AIR = 4;
PML = 5;
SIGMAINF = 6;
PEC = 7;
