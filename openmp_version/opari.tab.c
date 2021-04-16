#include "pomp_lib.h"


extern struct ompregdescr omp_rd_16;
extern struct ompregdescr omp_rd_17;
extern struct ompregdescr omp_rd_18;
extern struct ompregdescr omp_rd_19;
extern struct ompregdescr omp_rd_20;

int POMP_MAX_ID = 21;

struct ompregdescr* pomp_rd_table[21] = {
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  &omp_rd_16,
  &omp_rd_17,
  &omp_rd_18,
  &omp_rd_19,
  &omp_rd_20,
};
