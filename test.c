#include <string.h>

// tests

#define FFT_IMPLEMENTATION
#include "fft.h"

// test helpers
int failed = 0;
#define TEST(name) void name()
#define RUN_TEST(name) printf("\n\033[1m%s\n\033[0m", #name); name()
#define ASSERT(expr) if (!(expr)) { \
    failed = 1; \
    printf("\033[0;31mFAIL: %s\n\033[0m", #expr); \
} else { \
    printf("\033[0;32mPASS: %s\n\033[0m", #expr); \
}
#define ASSERT_STR_EQ(str1, str2) ASSERT(strcmp(str1, str2) == 0)
#define ASSERT_F_EQ(f1, f2) if ( \
    ((f1-f2)>-1.f*FFT_EPSILON) && \
    ((f1-f2)<FFT_EPSILON) \
    ) { \
    printf("\033[0;32mPASS: %s(%f)==%s(%f)\n\033[0m", #f1, f1, #f2, f2); \
} else { \
    printf("\033[0;31mFAIL: %s(%f)==%s(%f)\n\033[0m", #f1, f1, #f2, f2); \
}

int main() {
    return 0;
}
