#define Tb CONCAT2(Tb,nvec)
#define Y(arg) CONCAT2(arg,nvec)
#include "sharp_core_inc.c"
#if (MAXJOB_SPECIAL<6)
#include "sharp_core_inc3.c"
#endif

#if (MAXJOB_SPECIAL>=1)
#define njobs 1
#define Z(arg) CONCAT3(arg,nvec,njobs)
#include "sharp_core_inc2.c"
#undef Z
#undef njobs
#endif

#if (MAXJOB_SPECIAL>=2)
#define njobs 2
#define Z(arg) CONCAT3(arg,nvec,njobs)
#include "sharp_core_inc2.c"
#undef Z
#undef njobs
#endif

#if (MAXJOB_SPECIAL>=3)
#define njobs 3
#define Z(arg) CONCAT3(arg,nvec,njobs)
#include "sharp_core_inc2.c"
#undef Z
#undef njobs
#endif

#if (MAXJOB_SPECIAL>=4)
#define njobs 4
#define Z(arg) CONCAT3(arg,nvec,njobs)
#include "sharp_core_inc2.c"
#undef Z
#undef njobs
#endif

#if (MAXJOB_SPECIAL>=5)
#define njobs 5
#define Z(arg) CONCAT3(arg,nvec,njobs)
#include "sharp_core_inc2.c"
#undef Z
#undef njobs
#endif

#if (MAXJOB_SPECIAL>=6)
#define njobs 6
#define Z(arg) CONCAT3(arg,nvec,njobs)
#include "sharp_core_inc2.c"
#undef Z
#undef njobs
#endif

#undef Y
#undef Tb
