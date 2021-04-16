#include "pomp_lib.h"

struct ompregdescr omp_rd_16 = {
  "parallelfor", "", 0, "ballAlg.c", 111, 111, 121, 121
};

#define POMP_DLIST_00016 shared(omp_rd_16,omp_rd_17)

struct ompregdescr omp_rd_17 = {
  "critical", "find_a", 0, "ballAlg.c", 115, 115, 120, 120
};

struct ompregdescr omp_rd_18 = {
  "parallelfor", "", 0, "ballAlg.c", 126, 126, 136, 136
};

#define POMP_DLIST_00018 shared(omp_rd_18,omp_rd_19)

struct ompregdescr omp_rd_19 = {
  "critical", "find_b", 0, "ballAlg.c", 130, 130, 135, 135
};

struct ompregdescr omp_rd_20 = {
  "parallelfor", "", 0, "ballAlg.c", 152, 152, 171, 171
};

#define POMP_DLIST_00020 shared(omp_rd_20)

