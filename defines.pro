
colorpp = "Ivory";
colorro    = "LightGrey";

ppEM = "2Electromagnetic parameters/0";

mhz = 1.0e6;
speed_of_ligth = 299792458.0;
fopt = 1536.0;
fmin = 1530.0;
fmax = 1580.0;
nn = 15;

// DefineConstant[
//     freq = { fopt, Min fmin, Max fmax, Step (fmax - fmin) / nn,
//         Name StrCat[ppEM, "0Frequency [GHz]"], Loop 0, Highlight Str[colorpp],
//         Help Str["- Resonance frequency is 1.575 GHz'"]
//     },
//     lambda = {speed_of_ligth / (freq * mhz),
//         Name StrCat[ppEM, "1Wavelength [m]"], ReadOnly 1, Highlight Str[colorro]
//     },
//     k0 = {2.0 * Pi / lambda,
//         Name StrCat[ppEM, "2Wave number"], ReadOnly 1, Highlight Str[colorro]
//     }
// ];

freq = GetNumber["Model/Frequency"];
freq = freq * mhz;
lambda = speed_of_ligth / freq;
k0 = 2.0 * Pi / lambda;

SetNumber["Model/Lambda", lambda];
// Printf ["ololo: %g", GetNumber["Model/0Frequency"]];
