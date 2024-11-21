
#define FFT_IMPLEMENTATION
#include "fft.h"
#include <math.h>

#define N 10
#define EPSILON 0.00001f

void print_complex(float complex f) {
    printf("%.5f + %.5fi\n", creal(f), cimag(f));
}

void round_float(float* f) {
    if (*f < EPSILON && *f > -1.f*EPSILON) {
        *f = 0.f;
    }
}

void round_complex(float complex* f) {
    float re = creal(*f);
    float im = cimag(*f);

    round_float(&re);
    round_float(&im);

    *f = re + im*I;
}

int main()
{
    fft_Vec_cf a = fft_vec_alloc(N);
    fft_Vec_cf b = fft_vec_alloc(N);
#if 1
    for (size_t i=0; i<a.size; ++i) {
        fft_vec_set(a, i, (float)i);
        print_complex(fft_vec_get(a, i));
    }
#else
    fft_vec_set(a, 0, 1.f);
    fft_vec_set(a, 1, 2.f - 1.f*I);
    fft_vec_set(a, 2, -1.f*I);
    fft_vec_set(a, 3, -1.f + 2.f*I);
    for (size_t i=0; i<a.size; ++i) {
        round_complex(fft_vec_at(a, i));
        print_complex(fft_vec_get(a, i));
    }
#endif
    printf("----------------------\n");
    for (size_t i=0; i<b.size; ++i) {
        float complex res = 0.f;
        for (size_t j=0; j<a.size; ++j) {
            float complex xn = fft_vec_get(a, j);
            res += xn * cexp((-2.f*I*M_PI*(float)j*(float)i)/(float)b.size);
        }
        fft_vec_set(b, i, res);
    }
    for (size_t i=0; i<b.size; ++i) {
        print_complex(fft_vec_get(b, i));
    }

    printf("----------------------\n");
    for (size_t i=0; i<a.size; ++i) {
        float complex res = 0.f;
        for (size_t j=0; j<b.size; ++j) {
            float complex xn = fft_vec_get(b, j);
            res += xn * cexp((2.f*I*M_PI*(float)j*(float)i)/(float)a.size);
        }
        res /= a.size;
        fft_vec_set(a, i, res);
    }
    for (size_t i=0; i<a.size; ++i) {
        print_complex(fft_vec_get(a, i));
    }

    return 0;
}
