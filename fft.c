
#define FFT_IMPLEMENTATION
#include "fft.h"
#include <math.h>
#include <complex.h>

#define EPSILON 0.00001f
#define DEBUG 0
#define INFO if (DEBUG>0) printf


float complex omega(size_t n, size_t k)
{
    return cexp((-2.f*I*M_PI*k)/n);
}

void butterfly(fft_Vec_cf* series, fft_Vec_cf* dft_series)
{
    FFT_ASSERT(series->size == dft_series->size);
    if (series->size == 1) {
        fft_vec_set(dft_series, 0, fft_vec_get(series, 0));
        return;
    }
    size_t series_count = series->size;
    FFT_ASSERT(series_count % 2 == 0);

    // split the series into even and off
    fft_Vec_cf evn = fft_vec_alloc(series_count/2);
    fft_Vec_cf odd = fft_vec_alloc(series_count/2);
    for (size_t k = 0; k < series_count/2; ++k) {
        fft_vec_set(&evn, k, fft_vec_get(series, k*2));
        fft_vec_set(&odd, k, fft_vec_get(series, (k*2)+1));
    }

    // get dft of each part
    fft_Vec_cf evn_dft = fft_vec_alloc(series_count/2);
    fft_Vec_cf odd_dft = fft_vec_alloc(series_count/2);

    butterfly(&evn, &evn_dft);
    butterfly(&odd, &odd_dft);

    fft_vec_free(&evn);
    fft_vec_free(&odd);

    // butterfly them together
    for (size_t k = 0; k < series_count/2; ++k) {
        float complex w = omega(series_count, k);
        float complex x0 = fft_vec_get(&evn_dft, k);
        float complex x1 = fft_vec_get(&odd_dft, k);

        float complex y0 = x0 + x1*w;
        float complex y1 = x0 - x1*w;

        fft_vec_set(dft_series, k, y0);
        fft_vec_set(dft_series, k+series_count/2, y1);
    }

    fft_vec_free(&evn_dft);
    fft_vec_free(&odd_dft);
}

void print_complex(float complex f) {
    INFO("%.5f + %.5fi\n", creal(f), cimag(f));
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

#define N 1<<14
#define VAL_MAX 32

#define PERFORM_FFT

int main()
{
    fft_Vec_cf a = fft_vec_alloc(N);
    fft_Vec_cf b = fft_vec_alloc(N);
    fft_Vec_cf c = fft_vec_alloc(N);
    INFO("----------------------\nORIGINAL:\n");
#if 0
    for (size_t i=0; i<a.size; ++i) {
        fft_vec_set(&a, i, (float)(i % VAL_MAX));
        print_complex(fft_vec_get(&a, i));
    }
#else
    fft_vec_set(&a, 0, 1.f);
    fft_vec_set(&a, 1, 2.f - 1.f*I);
    fft_vec_set(&a, 2, -1.f*I);
    fft_vec_set(&a, 3, -1.f + 2.f*I);
    for (size_t i=0; i<a.size; ++i) {
        round_complex(fft_vec_at(&a, i));
        print_complex(fft_vec_get(&a, i));
    }
#endif

#ifdef PERFORM_FFT
    INFO("----------------------\nFFT:\n");
    butterfly(&a, &c);
    for (size_t i=0; i<c.size; ++i) {
        print_complex(fft_vec_get(&c, i));
    }
#else
    INFO("----------------------\nDFT:\n");
    for (size_t i=0; i<b.size; ++i) {
        float complex res = 0.f;
        for (size_t j=0; j<a.size; ++j) {
            float complex xn = fft_vec_get(&a, j);
            res += xn * cexp((-2.f*I*M_PI*(float)j*(float)i)/(float)b.size);
        }
        fft_vec_set(&b, i, res);
    }
    for (size_t i=0; i<b.size; ++i) {
        print_complex(fft_vec_get(&b, i));
    }
#endif

#if 0
    INFO("----------------------\nINVERSE DFT:\n");
    for (size_t i=0; i<a.size; ++i) {
        float complex res = 0.f;
        for (size_t j=0; j<b.size; ++j) {
            float complex xn = fft_vec_get(&b, j);
            res += xn * cexp((2.f*I*M_PI*(float)j*(float)i)/(float)a.size);
        }
        res /= a.size;
        fft_vec_set(&a, i, res);
    }
    for (size_t i=0; i<a.size; ++i) {
        print_complex(fft_vec_get(&a, i));
    }
#endif



    fft_vec_free(&a);
    fft_vec_free(&b);
    fft_vec_free(&c);
    
    return 0;
}
