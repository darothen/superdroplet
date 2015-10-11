
#ifndef SD_CPP_CONSTANTS_HPP
#define SD_CPP_CONSTANTS_HPP

namespace constants {

    enum resolution {
        lo = 2048, // 2^11
        med_lo = 8192, // 2^13
        med_hi = 32768, // 2^15
        hi = 131072, // 2^17
    };

    const double RHO_WATER = 1000.; // kg/m^3
}

#endif //SD_CPP_CONSTANTS_HPP
