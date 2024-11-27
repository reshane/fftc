
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

#define FFT_EPSILON 0.00001f

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

// [[ ALGORITHMS ]]

void dft(fft_Vec_cf* series, fft_Vec_cf* dft_series);
void dft_inverse(fft_Vec_cf* dft_series, fft_Vec_cf* series);
void fft(fft_Vec_cf* series, fft_Vec_cf* dft_series);
void fft_inverse(fft_Vec_cf* dft_series, fft_Vec_cf* series);

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

void fft_inverse(fft_Vec_cf* dft_series, fft_Vec_cf* series)
{
    FFT_ASSERT(dft_series->size == series->size);
    if (dft_series->size == 1) {
        fft_vec_set(series, 0, fft_vec_get(dft_series, 0));
        return;
    }
    size_t series_count = dft_series->size;
    FFT_ASSERT(series_count % 2 == 0);

    // split the series into even and off
    fft_Vec_cf evn = fft_vec_alloc(series_count/2);
    fft_Vec_cf odd = fft_vec_alloc(series_count/2);
    for (size_t k = 0; k < series_count/2; ++k) {
        fft_vec_set(&evn, k, fft_vec_get(dft_series, k*2));
        fft_vec_set(&odd, k, fft_vec_get(dft_series, (k*2)+1));
    }

    // get dft of each part
    fft_Vec_cf evn_dft = fft_vec_alloc(series_count/2);
    fft_Vec_cf odd_dft = fft_vec_alloc(series_count/2);

    fft_inverse(&evn, &evn_dft);
    fft_inverse(&odd, &odd_dft);

    fft_vec_free(&evn);
    fft_vec_free(&odd);

    // butterfly them together
    for (size_t k = 0; k < series_count/2; ++k) {
        float complex w = cexp((2.f*I*M_PI*(int)k)/series_count);
        float complex x0 = fft_vec_get(&evn_dft, k);
        float complex x1 = fft_vec_get(&odd_dft, k);

        float complex y0 = 0.5f*(x0 + x1*w);
        float complex y1 = 0.5f*(x0 - x1*w);

        fft_vec_set(series, k, y0);
        fft_vec_set(series, k+series_count/2, y1);
    }

    fft_vec_free(&evn_dft);
    fft_vec_free(&odd_dft);
}

void fft_matrix(fft_Matrix_cf* mat)
{
    // TODO: this should be its own standalone function
    // and should not rely on the fft of Vectors
    // should be able to apply fft to a row and col
    // also should accept a target matrix of the same dims
    size_t img_rows = mat->rows;
    size_t img_cols = mat->cols;
    fft_Vec_cf rows[img_rows];
    fft_Vec_cf dft_rows[img_rows];
    fft_Vec_cf cols[img_cols];
    fft_Vec_cf dft_cols[img_cols];

    // put the matrix rows in the rows array of vectors
    for (size_t k = 0; k < img_rows; ++k) {
        rows[k] = fft_vec_alloc(img_cols);
        for (size_t j = 0; j < img_cols; ++j) {
            fft_vec_set(&rows[k], j, fft_mat_get(mat, k, j));
        }
    }

    // perform an fft on each row
    for (size_t k = 0; k < img_rows; ++k) {
        dft_rows[k] = fft_vec_alloc(img_cols);
        fft(&rows[k], &dft_rows[k]);
        fft_vec_free(&rows[k]);
    }

    for (size_t k = 0; k < mat->rows; ++k) {
        for (size_t j = 0; j < mat->cols; ++j) {
            // float complex f = fft_vec_get(&dft_rows[k], j);
            // printf("%f+%fi ", creal(f), cimag(f));
        }
        // printf("\n");
    }

    // copy the rows over to the cols array of vectors
    for (size_t j = 0; j < img_cols; j++) {
        cols[j] = fft_vec_alloc(img_rows);
        for (size_t k = 0; k < img_rows; ++k) {
            fft_vec_set(&cols[j], k, fft_vec_get(&dft_rows[k], j));
        }
    }

    // free all the dft_rows vectors
    for (size_t k = 0; k < img_rows; ++k) {
        fft_vec_free(&dft_rows[k]);
    }

    // perform an fft on each column 
    for (size_t k = 0; k < img_cols; ++k) {
        dft_cols[k] = fft_vec_alloc(img_rows);
        fft(&cols[k], &dft_cols[k]);
        fft_vec_free(&cols[k]);
    }

    // copy the cols back into the matrix
    for (size_t k = 0; k < img_cols; ++k) {
        for (size_t j = 0; j < img_rows; ++j) {
            fft_mat_set(mat, j, k, fft_vec_get(&dft_cols[k], j));
        }
        fft_vec_free(&dft_cols[k]);
    }
}

void fft_matrix_inverse(fft_Matrix_cf* mat)
{
    size_t img_rows = mat->rows;
    size_t img_cols = mat->cols;
    fft_Vec_cf rows[img_rows];
    fft_Vec_cf dft_rows[img_rows];
    fft_Vec_cf cols[img_cols];
    fft_Vec_cf dft_cols[img_cols];

    // put the matrix cols in the cols array of vectors
    for (size_t k = 0; k < img_cols; ++k) {
        cols[k] = fft_vec_alloc(img_rows);
        for (size_t j = 0; j < img_rows; ++j) {
            fft_vec_set(&cols[k], j, fft_mat_get(mat, k, j));
        }
    }

    // perform an inverse fft on each col
    for (size_t k = 0; k < img_cols; ++k) {
        dft_cols[k] = fft_vec_alloc(img_rows);
        fft_inverse(&cols[k], &dft_cols[k]);
        fft_vec_free(&cols[k]);
    }

    // copy the cols over to the rows array of vectors
    for (size_t j = 0; j < img_rows; j++) {
        rows[j] = fft_vec_alloc(img_cols);
        for (size_t k = 0; k < img_cols; ++k) {
            fft_vec_set(&rows[j], k, fft_vec_get(&dft_cols[k], j));
        }
    }

    // free all the dft_cols vectors
    for (size_t k = 0; k < img_cols; ++k) {
        fft_vec_free(&dft_cols[k]);
    }

    // perform an fft on each row
    for (size_t k = 0; k < img_rows; ++k) {
        dft_rows[k] = fft_vec_alloc(img_cols);
        fft_inverse(&rows[k], &dft_rows[k]);
        fft_vec_free(&rows[k]);
    }

    // copy the rows back into the matrix
    for (size_t k = 0; k < img_rows; ++k) {
        for (size_t j = 0; j < img_cols; ++j) {
            fft_mat_set(mat, j, k, fft_vec_get(&dft_rows[k], j));
        }
        fft_vec_free(&dft_rows[k]);
    }
}

#endif // FFT_IMPLEMENTATION
