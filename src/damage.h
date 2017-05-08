/* headers */
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

/* constants */
#define PI      3.1415926535897932384626433832795
#define Cy      1740.0                      /* yield modulus, in MPa */
#define Dc      0.3                         /* damage critical value */
#define h       0.2                         /* closure parameter */
#define S       0.5                         /* damage parameter */
#define s       0.5                         /* damage parameter */
#define sig_f   180.0                       /* yield limit, in MPa */
#define sig_u   577.0                       /* endurance limit, in MPa */
#define eps_D   0.08                        /* ultimate strain */
#define wd      ((sig_u-sig_f)*eps_D)       /* damage initiation threshold */

