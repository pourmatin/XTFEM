/* headers */
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

/* constants */
#define PI      3.1415926535897932384626433832795
#define E       1.97e5                      /* Young's modulus, in MPa */
#define nu      0.3                         /* Poisson's ratio */
#define G       (E/(2.0*(1.0+nu)))          /* Shear modulus, in MPa */
#define K       (E/(3.0*(1.0-2.0*nu)))      /* bulk modulus, in MPa */
#define Cy      1740.0                      /* yield modulus, in MPa */
#define Dc      0.3                         /* damage critical value */
#define h       0.2                         /* closure parameter */
#define S       0.5                         /* damage parameter */
#define s       0.5                         /* damage parameter */
#define sig_f   180.0                       /* yield limit, in MPa */
#define sig_u   577.0                       /* endurance limit, in MPa */
#define eps_D   0.08                        /* ultimate strain */
#define a       ((1.0+nu)/(3.0*(1.0-nu)))   /* localization parameter */
#define b       (2.0*(4.0-5.0*nu)/(15.0*(1.0-nu)))  /* localization parameter */
#define wd      ((sig_u-sig_f)*eps_D)       /* damage initiation threshold */
#define COEF    (E/((1.0+nu)*(1.0-2.0*nu))) /* coefficient for stiffmat */

