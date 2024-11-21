
#ifndef FFT_H_
#define FFT_H_

// == INCLUDES ==

#include <stdio.h>
#include <complex.h>

#ifndef FFT_NO_UNI_STD
#include <unistd.h>
#endif // FFT_NO_UNI_STD

#ifndef FFT_MALLOC
#include <stdlib.h>
#define FFT_ALLOC malloc
#define FFT_FREE free
#endif // FFT_MALLOC

#ifndef FFT_ASSERT
#include <assert.h>
#define FFT_ASSERT assert
#endif // FFT_ASSERT


// == DEFINITIONS ==

// -- Vector --

// COMPLEX FLOAT

typedef struct {
    size_t size;
    float complex* data;
} fft_Vec_cf;
fft_Vec_cf fft_vec_alloc(size_t cap);
void fft_vec_free(fft_Vec_cf v);
float complex* fft_vec_at(fft_Vec_cf v, size_t i);
float complex fft_vec_get(fft_Vec_cf v, size_t i);
void fft_vec_set(fft_Vec_cf v, size_t i, float complex e);

#endif // FFT_H_

#ifdef FFT_IMPLEMENTATION

// == IMPLEMENTATION ==

// -- Vector --

// COMPLEX FLOAT

fft_Vec_cf fft_vec_alloc(size_t cap)
{
    fft_Vec_cf v;

    v.size = cap;
    v.data = FFT_ALLOC(v.size*sizeof(*v.data));
    FFT_ASSERT(v.data);

    return v;
}

void fft_vec_free(fft_Vec_cf v)
{
    FFT_FREE(v.data);
}

float complex* fft_vec_at(fft_Vec_cf v, size_t i)
{
    FFT_ASSERT(i < v.size);
    return &v.data[i];
}

float complex fft_vec_get(fft_Vec_cf v, size_t i)
{
    return *fft_vec_at(v, i);
}

void fft_vec_set(fft_Vec_cf v, size_t i, float complex e)
{
    *fft_vec_at(v, i) = e;
}

#endif // FFT_IMPLEMENTATION
