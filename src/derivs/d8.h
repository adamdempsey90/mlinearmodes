#define DACC  8
#define D1SIZE  9
#define D2SIZE  10

static const double d1cnum[D1SIZE][D1SIZE] = {{-761, 8, -14, 56, -35, 56, -14, 8, -1}, {-1, -223, 7, -7, 35, -7,
  7, -1, 1}, {1, -2, -19, 2, -5, 2, -1, 2, -1}, {-1, 1, -1, -9, 5, -1,
   1, -1, 1}, {1, -4, 1, -4, 0, 4, -1, 4, -1}, {-1, 1, -1, 1, -5, 9,
  1, -1, 1}, {1, -2, 1, -2, 5, -2, 19, 2, -1}, {-1, 1, -7, 7, -35,
  7, -7, 223, 1}, {1, -8, 14, -56, 35, -56, 14, -8, 761}};
static const double d1cden[D1SIZE][D1SIZE] = {{280, 1, 1, 3, 2, 5, 3, 7, 8}, {8, 140, 2, 2, 12, 4, 10, 6, 56}, {56,
   7, 20, 1, 4, 3, 4, 35, 168}, {168, 14, 2, 20, 4, 2, 6, 28,
  280}, {280, 105, 5, 5, 1, 5, 5, 105, 280}, {280, 28, 6, 2, 4, 20, 2,
   14, 168}, {168, 35, 4, 3, 4, 1, 20, 7, 56}, {56, 6, 10, 4, 12, 2,
  2, 140, 8}, {8, 7, 3, 5, 2, 3, 1, 1, 280}};
static const int d1off[D1SIZE][D1SIZE] = {{0, 1, 2, 3, 4, 5, 6, 7, 8}, {-1, 0, 1, 2, 3, 4, 5, 6, 7}, {-2, -1,
  0, 1, 2, 3, 4, 5, 6}, {-3, -2, -1, 0, 1, 2, 3, 4,
  5}, {-4, -3, -2, -1, 0, 1, 2, 3, 4}, {-5, -4, -3, -2, -1, 0, 1, 2,
  3}, {-6, -5, -4, -3, -2, -1, 0, 1, 2}, {-7, -6, -5, -4, -3, -2, -1,
  0, 1}, {-8, -7, -6, -5, -4, -3, -2, -1, 0}};

static const double d2cnum[D1SIZE][D2SIZE] = {{6515, -4609, 5869, -6289, 6499, -265, 6709, -967, 3407, -761}, {761,
   61, -201, 341, -1163, 411, -17, 1303, -9, 223}, {-223,
  293, -395, -13, 83, -319, 59, -5, 389, -19}, {19, -67, 97, -89, 23,
  7, -17, 11, -1, 1}, {-1, 8, -1, 8, -205, 8, -1, 8, -1, 0}, {1, -1,
  11, -17, 7, 23, -89, 97, -67, 19}, {-19, 389, -5, 59, -319,
  83, -13, -395, 293, -223}, {223, -9, 1303, -17, 411, -1163,
  341, -201, 61, 761}, {-761, 3407, -967, 6709, -265, 6499, -6289,
  5869, -4609, 6515}};
static const double d2cden[D1SIZE][D2SIZE] = {{1008, 140, 70, 45, 40, 2, 90, 35, 560, 1260}, {1260, 144, 35, 30,
  90, 40, 3, 630, 20, 5040}, {5040, 280, 252, 30, 40, 180, 60, 14,
  5040, 2520}, {2520, 560, 70, 36, 20, 40, 90, 140, 56, 560}, {560,
  315, 5, 5, 72, 5, 5, 315, 560, 1}, {560, 56, 140, 90, 40, 20, 36,
  70, 560, 2520}, {2520, 5040, 14, 60, 180, 40, 30, 252, 280,
  5040}, {5040, 20, 630, 3, 40, 90, 30, 35, 144, 1260}, {1260, 560,
  35, 90, 2, 40, 45, 70, 140, 1008}};
static const int d2off[D1SIZE][D2SIZE] = {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, {-1, 0, 1, 2, 3, 4, 5, 6, 7,
  8}, {-2, -1, 0, 1, 2, 3, 4, 5, 6, 7}, {-3, -2, -1, 0, 1, 2, 3, 4, 5,
   6}, {-4, -3, -2, -1, 0, 1, 2, 3, 4, 5}, {-6, -5, -4, -3, -2, -1, 0,
   1, 2, 3}, {-7, -6, -5, -4, -3, -2, -1, 0, 1,
  2}, {-8, -7, -6, -5, -4, -3, -2, -1, 0,
  1}, {-9, -8, -7, -6, -5, -4, -3, -2, -1, 0}};
