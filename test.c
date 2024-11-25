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
    failed = 1; \
}
#define ASSERT_CF_EQ(cf1, cf2) if ( \
    ((creal(cf1)-creal(cf2))>-1.f*FFT_EPSILON) && ((creal(cf1)-creal(cf2))<FFT_EPSILON) \
    && \
    ((cimag(cf1)-cimag(cf2))>-1.f*FFT_EPSILON) && ((cimag(cf1)-cimag(cf2))<FFT_EPSILON) \
    ) { \
    printf("\033[0;32mPASS: %s(%0.4f+%0.4fi)==%s(%0.4f+%0.4fi)\n\033[0m", \
           #cf1, creal(cf1), cimag(cf1), #cf2, creal(cf2), cimag(cf2)); \
} else { \
    printf("\033[0;31mFAIL: %s(%0.9f+%0.9fi)==%s(%0.9f+%0.9fi)\n\033[0m", \
           #cf1, creal(cf1), cimag(cf1), #cf2, creal(cf2), cimag(cf2)); \
    failed =  1; \
}
#define ASSERT_CF_EQ_S(cf1, cf2) if (!( \
    ((creal(cf1)-creal(cf2))>-1.f*FFT_EPSILON) && ((creal(cf1)-creal(cf2))<FFT_EPSILON) \
    && \
    ((cimag(cf1)-cimag(cf2))>-1.f*FFT_EPSILON) && ((cimag(cf1)-cimag(cf2))<FFT_EPSILON) \
    )) { \
    printf("\033[0;31mFAIL: %s(%0.9f+%0.9fi)==%s(%0.9f+%0.9fi)\n\033[0m", \
           #cf1, creal(cf1), cimag(cf1), #cf2, creal(cf2), cimag(cf2)); \
    failed =  1; \
}

TEST(test_dft)
{
    size_t n = 8;
    fft_Vec_cf a = fft_vec_alloc(n);
    fft_Vec_cf b = fft_vec_alloc(n);
    fft_Vec_cf c = fft_vec_alloc(n);

    fft_vec_set(&a, 0,  1.f +  0.f * I);
    fft_vec_set(&a, 1,  2.f + -1.f * I);
    fft_vec_set(&a, 2,  0.f + -1.f * I);
    fft_vec_set(&a, 3, -1.f +  2.f * I);
    fft_vec_set(&a, 4,  1.f +  0.f * I);
    fft_vec_set(&a, 5,  2.f + -1.f * I);
    fft_vec_set(&a, 6,  0.f + -1.f * I);
    fft_vec_set(&a, 7, -1.f +  2.f * I);


    fft_vec_set(&c, 0,  4.f +  0.f * I);
    fft_vec_set(&c, 1,  0.f +  0.f * I);
    fft_vec_set(&c, 2, -4.f + -4.f * I);
    fft_vec_set(&c, 3,  0.f +  0.f * I);
    fft_vec_set(&c, 4,  0.f + -4.f * I);
    fft_vec_set(&c, 5,  0.f +  0.f * I);
    fft_vec_set(&c, 6,  8.f +  8.f * I);
    fft_vec_set(&c, 7,  0.f +  0.f * I);

    dft(&a, &b);
    fft_vec_free(&a);

    for (size_t k = 0; k < n; ++k) {
        ASSERT_CF_EQ(fft_vec_get(&b, k), fft_vec_get(&c, k))
    }

    fft_vec_free(&b);
    fft_vec_free(&c);
}

TEST(test_dft_inverse)
{
    size_t n = 8;
    fft_Vec_cf a = fft_vec_alloc(n);
    fft_Vec_cf b = fft_vec_alloc(n);
    fft_Vec_cf c = fft_vec_alloc(n);

    fft_vec_set(&a, 0,  1.f +  0.f * I);
    fft_vec_set(&a, 1,  2.f + -1.f * I);
    fft_vec_set(&a, 2,  0.f + -1.f * I);
    fft_vec_set(&a, 3, -1.f +  2.f * I);
    fft_vec_set(&a, 4,  1.f +  0.f * I);
    fft_vec_set(&a, 5,  2.f + -1.f * I);
    fft_vec_set(&a, 6,  0.f + -1.f * I);
    fft_vec_set(&a, 7, -1.f +  2.f * I);


    fft_vec_set(&c, 0,  4.f +  0.f * I);
    fft_vec_set(&c, 1,  0.f +  0.f * I);
    fft_vec_set(&c, 2, -4.f + -4.f * I);
    fft_vec_set(&c, 3,  0.f +  0.f * I);
    fft_vec_set(&c, 4,  0.f + -4.f * I);
    fft_vec_set(&c, 5,  0.f +  0.f * I);
    fft_vec_set(&c, 6,  8.f +  8.f * I);
    fft_vec_set(&c, 7,  0.f +  0.f * I);

    dft_inverse(&c, &b);
    fft_vec_free(&c);

    for (size_t k = 0; k < n; ++k) {
        ASSERT_CF_EQ(fft_vec_get(&a, k), fft_vec_get(&b, k))
    }
    fft_vec_free(&a);
    fft_vec_free(&b);
}

TEST(test_fft)
{
    size_t n = 8;
    fft_Vec_cf a = fft_vec_alloc(n);
    fft_Vec_cf b = fft_vec_alloc(n);
    fft_Vec_cf c = fft_vec_alloc(n);

    fft_vec_set(&a, 0,  1.f +  0.f * I);
    fft_vec_set(&a, 1,  2.f + -1.f * I);
    fft_vec_set(&a, 2,  0.f + -1.f * I);
    fft_vec_set(&a, 3, -1.f +  2.f * I);
    fft_vec_set(&a, 4,  1.f +  0.f * I);
    fft_vec_set(&a, 5,  2.f + -1.f * I);
    fft_vec_set(&a, 6,  0.f + -1.f * I);
    fft_vec_set(&a, 7, -1.f +  2.f * I);


    fft_vec_set(&c, 0,  4.f +  0.f * I);
    fft_vec_set(&c, 1,  0.f +  0.f * I);
    fft_vec_set(&c, 2, -4.f + -4.f * I);
    fft_vec_set(&c, 3,  0.f +  0.f * I);
    fft_vec_set(&c, 4,  0.f + -4.f * I);
    fft_vec_set(&c, 5,  0.f +  0.f * I);
    fft_vec_set(&c, 6,  8.f +  8.f * I);
    fft_vec_set(&c, 7,  0.f +  0.f * I);

    fft(&a, &b);
    fft_vec_free(&a);

    for (size_t k = 0; k < n; ++k) {
        ASSERT_CF_EQ(fft_vec_get(&b, k), fft_vec_get(&c, k))
    }
    fft_vec_free(&b);
    fft_vec_free(&c);
}

TEST(test_fft_inverse)
{
    size_t n = 8;
    fft_Vec_cf a = fft_vec_alloc(n);
    fft_Vec_cf b = fft_vec_alloc(n);
    fft_Vec_cf c = fft_vec_alloc(n);

    fft_vec_set(&a, 0,  1.f +  0.f * I);
    fft_vec_set(&a, 1,  2.f + -1.f * I);
    fft_vec_set(&a, 2,  0.f + -1.f * I);
    fft_vec_set(&a, 3, -1.f +  2.f * I);
    fft_vec_set(&a, 4,  1.f +  0.f * I);
    fft_vec_set(&a, 5,  2.f + -1.f * I);
    fft_vec_set(&a, 6,  0.f + -1.f * I);
    fft_vec_set(&a, 7, -1.f +  2.f * I);


    fft_vec_set(&c, 0,  4.f +  0.f * I);
    fft_vec_set(&c, 1,  0.f +  0.f * I);
    fft_vec_set(&c, 2, -4.f + -4.f * I);
    fft_vec_set(&c, 3,  0.f +  0.f * I);
    fft_vec_set(&c, 4,  0.f + -4.f * I);
    fft_vec_set(&c, 5,  0.f +  0.f * I);
    fft_vec_set(&c, 6,  8.f +  8.f * I);
    fft_vec_set(&c, 7,  0.f +  0.f * I);

    fft_inverse(&c, &b);
    fft_vec_free(&c);

    for (size_t k = 0; k < n; ++k) {
        ASSERT_CF_EQ(fft_vec_get(&a, k), fft_vec_get(&b, k))
    }
    fft_vec_free(&a);
    fft_vec_free(&b);
}

TEST(test_dft_polynomial_multiplication)
{
    size_t n = 4;
    // (x^2 + x)(x^2 + x) = (x^4 + 2x^3 + x^2 0x)
    fft_Vec_cf a = fft_vec_alloc(n);
    fft_Vec_cf b = fft_vec_alloc(n);
    fft_Vec_cf c = fft_vec_alloc(n);
    fft_Vec_cf d = fft_vec_alloc(n);

    fft_vec_set(&a, 0, 0.f);
    fft_vec_set(&a, 1, 0.f);
    fft_vec_set(&a, 2, 1.f);
    fft_vec_set(&a, 3, 1.f);

    fft_vec_set(&b, 0, 0.f);
    fft_vec_set(&b, 1, 0.f);
    fft_vec_set(&b, 2, 1.f);
    fft_vec_set(&b, 3, 1.f);

    dft(&a, &c);
    fft_vec_free(&a);
    dft(&b, &d);
    fft_vec_free(&b);

    fft_Vec_cf e = fft_vec_alloc(n);
    fft_Vec_cf f = fft_vec_alloc(n);
    for (size_t k = 0; k < n; ++k) {
        fft_vec_set(&e, k, fft_vec_get(&c, k) * fft_vec_get(&d, k));
    }
    fft_vec_free(&c);
    fft_vec_free(&d);

    dft_inverse(&e, &f);
    fft_vec_free(&e);

    float complex g[n];
    g[0] = 1.f;
    g[1] = 2.f;
    g[2] = 1.f;
    g[3] = 0.f;

    for (size_t k = 0; k < n; ++k) {
        ASSERT_CF_EQ(fft_vec_get(&f, k), g[k])
    }
    fft_vec_free(&f);
}

TEST(test_fft_polynomial_multiplication)
{
    size_t n = 4;
    // (x^2 + x)(x^2 + x) = (x^4 + 2x^3 + x^2 0x)
    fft_Vec_cf a = fft_vec_alloc(n);
    fft_Vec_cf b = fft_vec_alloc(n);
    fft_Vec_cf c = fft_vec_alloc(n);
    fft_Vec_cf d = fft_vec_alloc(n);

    fft_vec_set(&a, 0, 0.f);
    fft_vec_set(&a, 1, 0.f);
    fft_vec_set(&a, 2, 1.f);
    fft_vec_set(&a, 3, 1.f);

    fft_vec_set(&b, 0, 0.f);
    fft_vec_set(&b, 1, 0.f);
    fft_vec_set(&b, 2, 1.f);
    fft_vec_set(&b, 3, 1.f);

    fft(&a, &c);
    fft_vec_free(&a);
    fft(&b, &d);
    fft_vec_free(&b);

    fft_Vec_cf e = fft_vec_alloc(n);
    fft_Vec_cf f = fft_vec_alloc(n);
    for (size_t k = 0; k < n; ++k) {
        fft_vec_set(&e, k, fft_vec_get(&c, k) * fft_vec_get(&d, k));
    }
    fft_vec_free(&c);
    fft_vec_free(&d);

    fft_inverse(&e, &f);
    fft_vec_free(&e);

    float complex g[n];
    g[0] = 1.f;
    g[1] = 2.f;
    g[2] = 1.f;
    g[3] = 0.f;

    for (size_t k = 0; k < n; ++k) {
        ASSERT_CF_EQ(fft_vec_get(&f, k), g[k])
    }
    fft_vec_free(&f);
}

TEST(fft_and_dft)
{
    size_t n = 1<<5;
    size_t val_max = 32;
    fft_Vec_cf a = fft_vec_alloc(n);
    fft_Vec_cf b = fft_vec_alloc(n);
    fft_Vec_cf c = fft_vec_alloc(n);
    fft_Vec_cf d = fft_vec_alloc(n);
    fft_Vec_cf e = fft_vec_alloc(n);

    for (size_t i=0; i<a.size; ++i) {
        fft_vec_set(&a, i, (float)(i % val_max));
    }

    fft(&a, &b);
    fft_inverse(&b, &d);
    fft_vec_free(&b);
    dft(&a, &c);
    dft_inverse(&c, &e);
    fft_vec_free(&c);

    printf("------------------------------\n");
    printf("|            \033[1mFFT\033[0m             |\n");
    printf("------------------------------\n");
    for (size_t k = 0; k < n; ++k) {
        ASSERT_CF_EQ_S(fft_vec_get(&a, k), fft_vec_get(&d, k))
    }
    printf("------------------------------\n");
    printf("|            \033[1mDFT\033[0m             |\n");
    printf("------------------------------\n");
    for (size_t k = 0; k < n; ++k) {
        ASSERT_CF_EQ_S(fft_vec_get(&a, k), fft_vec_get(&e, k))
    }
    fft_vec_free(&a);
    fft_vec_free(&d);
    fft_vec_free(&e);
}

int main() {
    RUN_TEST(test_dft);
    RUN_TEST(test_dft_inverse);
    RUN_TEST(test_dft_polynomial_multiplication);
    RUN_TEST(test_fft);
    RUN_TEST(test_fft_inverse);
    RUN_TEST(test_fft_polynomial_multiplication);
    RUN_TEST(fft_and_dft);
    return failed;
}
