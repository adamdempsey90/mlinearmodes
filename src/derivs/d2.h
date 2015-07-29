#define DACC  2
#define D1SIZE  3
#define D2SIZE  4

static const double d1cnum[D1SIZE][D1SIZE] = {{-3, 2, -1}, {-1, 0, 1}, {1, -2, 3}};
static const double d1cden[D1SIZE][D1SIZE] = {{2, 1, 2}, {2, 1, 2}, {2, 1, 2}};
static const int d1off[D1SIZE][D1SIZE] = {{0, 1, 2}, {-1, 0, 1}, {-2, -1, 0}};

static const double d2cnum[D1SIZE][D2SIZE] = {{2, -5, 4, -1}, {1, -2, 1, 0}, {-1, 4, -5, 2}};
static const double d2cden[D1SIZE][D2SIZE] = {{1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}};
static const int d2off[D1SIZE][D2SIZE] = {{0, 1, 2, 3}, {-1, 0, 1, 2}, {-3, -2, -1, 0}};
