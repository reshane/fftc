
#ifndef FFT_H_
#define FFT_H_

// == INCLUDES ==

#include <stdio.h>
#include <math.h>
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

// [[ CONSTANTS ]]

#define FFT_EPSILON 0.0000001f

// [[ STRUCTURES ]]

// -- Vector --

typedef struct {
    size_t size;
    float complex* data;
} fft_Vec_cf;
fft_Vec_cf fft_vec_alloc(size_t cap);
void fft_vec_free(fft_Vec_cf* v);
// deprecated
float complex* fft_vec_at(fft_Vec_cf* v, size_t i);
float complex fft_vec_get(fft_Vec_cf* v, size_t i);
void fft_vec_set(fft_Vec_cf* v, size_t i, float complex e);

// -- Matrix --
typedef struct {
    size_t rows;
    size_t cols;
    float complex* data;
} fft_Matrix_cf;
fft_Matrix_cf fft_mat_alloc(size_t rows, size_t cols);
void fft_mat_free(fft_Matrix_cf* mat);
float complex fft_mat_get(fft_Matrix_cf* mat, size_t x, size_t y);
void fft_mat_set(fft_Matrix_cf* mat, size_t x, size_t y, float complex e);

#endif // FFT_H_

#ifdef FFT_IMPLEMENTATION

// == IMPLEMENTATION ==

// [[ STRUCTURES ]]

// -- Vector --

fft_Vec_cf fft_vec_alloc(size_t cap)
{
    fft_Vec_cf v;

    v.size = cap;
    v.data = FFT_ALLOC(v.size*sizeof(*v.data));
    FFT_ASSERT(v.data);

    return v;
}

void fft_vec_free(fft_Vec_cf* v)
{
    FFT_FREE(v->data);
}

float complex* fft_vec_at(fft_Vec_cf* v, size_t i)
{
    FFT_ASSERT(i < v->size);
    return &v->data[i];
}

float complex fft_vec_get(fft_Vec_cf* v, size_t i)
{
    return *fft_vec_at(v, i);
}

void fft_vec_set(fft_Vec_cf* v, size_t i, float complex e)
{
    FFT_ASSERT(i < v->size);
    *fft_vec_at(v, i) = e;
}

// -- Matrix --

fft_Matrix_cf fft_mat_alloc(size_t rows, size_t cols)
{
    fft_Matrix_cf mat;
    mat.rows = rows;
    mat.cols = cols;
    mat.data = (float complex*)malloc(sizeof(float complex)*rows*cols);
    FFT_ASSERT(mat.data);
    return mat;
}

void fft_mat_free(fft_Matrix_cf* mat)
{
    free(mat->data);
}

float complex fft_mat_get(fft_Matrix_cf* mat, size_t x, size_t y)
{
    size_t idx = y*mat->cols + x;
    FFT_ASSERT(idx<(mat->cols)*(mat->rows));
    return mat->data[idx];
}

void fft_mat_set(fft_Matrix_cf* mat, size_t x, size_t y, float complex e)
{
    size_t idx = y*mat->cols + x;
    FFT_ASSERT(idx<(mat->cols)*(mat->rows));
    mat->data[idx] = e;
}

// [[ ALGORITHMS ]]

// -- Discrete Fourier Transform -- 

void dft(fft_Vec_cf* series, fft_Vec_cf* dft_series)
{
    for (size_t i=0; i<dft_series->size; ++i) {
        float complex res = 0.f;
        for (size_t j=0; j<series->size; ++j) {
            float complex xn = fft_vec_get(series, j);
            res += xn * cexp((-2.f*I*M_PI*(float)j*(float)i)/(float)dft_series->size);
        }
        fft_vec_set(dft_series, i, res);
    }
}

void dft_inverse(fft_Vec_cf* dft_series, fft_Vec_cf* series)
{
    for (size_t i=0; i<series->size; ++i) {
        float complex res = 0.f;
        for (size_t j=0; j<dft_series->size; ++j) {
            float complex xn = fft_vec_get(dft_series, j);
            res += xn * cexp((2.f*I*M_PI*(float)j*(float)i)/(float)series->size);
        }
        res /= series->size;
        fft_vec_set(series, i, res);
    }
}

// -- Fast Fourier Transform --

float complex fft_omega(size_t n, size_t k)
{
    return cexp((-2.f*I*M_PI*k)/n);
}

void fft(fft_Vec_cf* series, fft_Vec_cf* dft_series)
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

    fft(&evn, &evn_dft);
    fft(&odd, &odd_dft);

    fft_vec_free(&evn);
    fft_vec_free(&odd);

    // butterfly them together
    for (size_t k = 0; k < series_count/2; ++k) {
        float complex w = fft_omega(series_count, k);
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

#endif // FFT_IMPLEMENTATION
