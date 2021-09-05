
mhz = 1.0e6;
speed_of_ligth = 299792458.0;
freq = GetNumber["Model/Frequency"];
freq = freq * mhz;
lambda = speed_of_ligth / freq;
k0 = 2.0 * Pi / lambda;
epr = GetNumber["Model/epr"];

SetNumber["Model/Lambda", lambda];
SetNumber["Model/WaveNumber", k0];
// Printf["ololo: %g", GetNumber["Model/Frequency"]];
Printf["ololo!"];
Printf[OnelabAction];
