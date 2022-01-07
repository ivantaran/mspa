
speed_of_ligth = 299792458.0;

fmin_mhz = 400.0;
fmax_mhz = 450.0;
mhz = 1.0e6;
freq_step_num = 10;
freq_step = (fmax_mhz - fmin_mhz) / freq_step_num;

DefineConstant[
    freq = {
        fmin_mhz + 35.0, 
        Min fmin_mhz, 
        Max fmax_mhz, 
        Step freq_step,
        Name "Model/0Frequency [MHz]"
    },
    lambda = { speed_of_ligth / freq / mhz,
        Name "Model/1Wavelength [m]", 
        ReadOnly 1
    },
    k0 = {
        2.0 * Pi / lambda,
        Name "Model/2Wave number", 
        ReadOnly 1
    }
] ;
// (atan(sqrt(2.0)*y/z)+(1.0-sign(z))*sign(y)*pi*0.5)*180.0/pi  
// freq = GetNumber["freq"];
// freq = GetNumber["Model/Frequency"];
// freq = 140.0;
freq = freq * mhz;
// lambda = speed_of_ligth / freq;
// k0 = 2.0 * Pi / lambda;
// epr = GetNumber["Model/epr"];
// epr = 1.5;

// SetNumber["Model/Lambda", lambda];
// SetNumber["Model/WaveNumber", k0];
// Printf["ololo: %g", GetNumber["Model/Frequency"]];
Printf["ololo!"];
Printf[OnelabAction];

SKINFEED = 1;
SKINCONDUCTOR = 2;
AIR = 4;
PML = 5;
SIGMAINF = 6;
PEC = 7;
