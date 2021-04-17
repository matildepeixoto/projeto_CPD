#include "pomp_lib.h"

struct ompregdescr omp_rd_121 = {
  "parallelfor", "", 0, "ballAlg.c", 109, 109, 119, 119
};

#define POMP_DLIST_00121 shared(omp_rd_121,omp_rd_122)

struct ompregdescr omp_rd_122 = {
  "critical", "find_a", 0, "ballAlg.c", 113, 113, 118, 118
};

struct ompregdescr omp_rd_123 = {
  "parallelfor", "", 0, "ballAlg.c", 124, 124, 134, 134
};

#define POMP_DLIST_00123 shared(omp_rd_123,omp_rd_124)

struct ompregdescr omp_rd_124 = {
  "critical", "find_b", 0, "ballAlg.c", 128, 128, 133, 133
};

struct ompregdescr omp_rd_125 = {
  "parallelfor", "", 0, "ballAlg.c", 150, 150, 169, 169
};

#define POMP_DLIST_00125 shared(omp_rd_125)

struct ompregdescr omp_rd_126 = {
  "parallelfor", "", 0, "ballAlg.c", 262, 262, 269, 269
};

#define POMP_DLIST_00126 shared(omp_rd_126,omp_rd_127)

struct ompregdescr omp_rd_127 = {
  "critical", "find_radius", 0, "ballAlg.c", 266, 266, 268, 268
};

struct ompregdescr omp_rd_128 = {
  "parallelsections", "", 2, "ballAlg.c", 320, 320, 326, 326
};

#define POMP_DLIST_00128 shared(omp_rd_128)

struct ompregdescr omp_rd_129 = {
  "critical", "", 0, "ballAlg.c", 329, 329, 330, 330
};

struct ompregdescr omp_rd_130 = {
  "parallel", "", 0, "ballAlg.c", 332, 332, 341, 341
};

#define POMP_DLIST_00130 shared(omp_rd_130,omp_rd_131)

struct ompregdescr omp_rd_131 = {
  "sections", "", 2, "ballAlg.c", 334, 334, 340, 340
};

struct ompregdescr omp_rd_132 = {
  "critical", "", 0, "ballAlg.c", 350, 350, 351, 351
};

