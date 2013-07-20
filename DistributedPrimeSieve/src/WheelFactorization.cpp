#include "WheelFactorization.h"

const WheelElement wheel30Elements[30] = {
/* spoke 1    */{ 0xFF, 1 },
/* spoke 7    */{ 0, 6 }, { 0xFF, 5 }, { 0xFF, 4 }, { 0xFF, 3 }, { 0xFF, 2 }, { 0xFF, 1 },
/* spoke 11   */{ 1, 4 }, { 0xFF, 3 }, { 0xFF, 2 }, { 0xFF, 1 },
/* spoke 13   */{ 2, 2 }, { 0xFF, 1 },
/* spoke 17   */{ 3, 4 }, { 0xFF, 3 }, { 0xFF, 2 }, { 0xFF, 1 },
/* spoke 19   */{ 4, 2 }, { 0xFF, 1 },
/* spoke 23   */{ 5, 4 }, { 0xFF, 3 }, { 0xFF, 2 }, { 0xFF, 1 },
/* spoke 29   */{ 6, 6 }, { 0xFF, 5 }, { 0xFF, 4 }, { 0xFF, 3 }, { 0xFF, 2 }, { 0xFF, 1 },
/* next spoke */{ 7, 2 } };

const WheelElement wheel210Elements[210] = {
/* spoke 1    */{  0xFF, 1 },
/* spoke 11   */{  0, 10 }, { 0xFF, 9  }, { 0xFF, 8  }, { 0xFF, 7  }, { 0xFF, 6  }, { 0xFF, 5  }, { 0xFF, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* spoke 13   */{  1, 2  }, { 0xFF, 1  },
/* spoke 17   */{  2, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* spoke 19   */{  3, 2  }, { 0xFF, 1  },
/* spoke 23   */{  4, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* spoke 29   */{  5, 6  }, { 0xFF, 5  }, { 0xFF, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* spoke 31   */{  6, 2  }, { 0xFF, 1  },
/* spoke 37   */{  7, 6  }, { 0xFF, 5  }, { 0xFF, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* spoke 41   */{  8, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* spoke 43   */{  9, 2  }, { 0xFF, 1  },
/* spoke 47   */{ 10, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* spoke 53   */{ 11, 6  }, { 0xFF, 5  }, { 0xFF, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* spoke 59   */{ 12, 6  }, { 0xFF, 5  }, { 0xFF, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* spoke 61   */{ 13, 2  }, { 0xFF, 1  },
/* spoke 67   */{ 14, 6  }, { 0xFF, 5  }, { 0xFF, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* spoke 71   */{ 15, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* spoke 73   */{ 16, 2  }, { 0xFF, 1  },
/* spoke 79   */{ 17, 6  }, { 0xFF, 5  }, { 0xFF, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* spoke 83   */{ 18, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* spoke 89   */{ 19, 6  }, { 0xFF, 5  }, { 0xFF, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* spoke 97   */{ 20, 8  }, { 0xFF, 7  }, { 0xFF, 6  }, { 0xFF, 5  }, { 0xFF, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* spoke 101  */{ 21, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* spoke 103  */{ 22, 2  }, { 0xFF, 1  },
/* spoke 107  */{ 23, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* spoke 109  */{ 24, 2  }, { 0xFF, 1  },
/* spoke 113  */{ 25, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* spoke 121  */{ 26, 8  }, { 0xFF, 7  }, { 0xFF, 6  }, { 0xFF, 5  }, { 0xFF, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* spoke 127  */{ 27, 6  }, { 0xFF, 5  }, { 0xFF, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* spoke 131  */{ 28, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* spoke 137  */{ 29, 6  }, { 0xFF, 5  }, { 0xFF, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* spoke 139  */{ 30, 2  }, { 0xFF, 1  },
/* spoke 143  */{ 31, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* spoke 149  */{ 32, 6  }, { 0xFF, 5  }, { 0xFF, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* spoke 151  */{ 33, 2  }, { 0xFF, 1  },
/* spoke 157  */{ 34, 6  }, { 0xFF, 5  }, { 0xFF, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* spoke 163  */{ 35, 6  }, { 0xFF, 5  }, { 0xFF, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* spoke 167  */{ 36, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* spoke 169  */{ 37, 2  }, { 0xFF, 1  },
/* spoke 173  */{ 38, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* spoke 179  */{ 39, 6  }, { 0xFF, 5  }, { 0xFF, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* spoke 181  */{ 40, 2  }, { 0xFF, 1  },
/* spoke 187  */{ 41, 6  }, { 0xFF, 5  }, { 0xFF, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* spoke 191  */{ 42, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* spoke 193  */{ 43, 2  }, { 0xFF, 1  },
/* spoke 197  */{ 44, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* spoke 199  */{ 45, 2  }, { 0xFF, 1  },
/* spoke 209  */{ 46, 10 }, { 0xFF, 9  }, { 0xFF, 8  }, { 0xFF, 7  }, { 0xFF, 6  }, { 0xFF, 5  }, { 0xFF, 4  }, { 0xFF, 3  }, { 0xFF, 2  }, { 0xFF, 1  },
/* next spoke */{ 47, 2  } };