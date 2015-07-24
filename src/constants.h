#if !defined (CONSTANTS_H)
#define CONSTANTS_H


namespace PIC {

/*********************
* Physical constants *
*********************/
class Constants
{
public:
    static double c()       { return 2.9979e+10; } // light velocity in vacuum (cm/sec)
    static double e()       { return 4.8032e-10; } // electron charge
    static double mp()      { return 1.6726e-24; } // proton mass
    static double pi()      { return 3.14159265358979323846; } // pi
    static double e_mp()    { return 2.8717e+14; } // proton charge/mass ratio
    static double c_4pi_e() { return 4.96679926e+18; } // c/(4*pi*e)
};

} // namespace PIC

#endif //CONSTANTS_H
