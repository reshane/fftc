
#define FFT_IMPLEMENTATION
#include "fft.h"
#include <stdint.h>
#include <errno.h>
// #include <strerror.h>

void fft_matrix(fft_Matrix_cf* mat)
{
    // TODO: this should be its own standalone function
    // and should not rely on the fft of Vectors
    // should be able to apply fft to a row and col
    // also should accept a target matrix of the same dims
    size_t img_rows = mat->rows;
    size_t img_cols = mat->cols;
    printf("rows: %ld, cols: %ld\n", img_rows, img_cols);
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
    printf("rows: %ld, cols: %ld\n", img_rows, img_cols);
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

void print_matrix(fft_Matrix_cf* mat)
{
    for (size_t k = 0; k < mat->rows; ++k) {
        for (size_t j = 0; j < mat->cols; ++j) {
            printf("%f+%fi ", creal(fft_mat_get(mat, k, j)), cimag(fft_mat_get(mat, k, j)));
        }
        printf("\n");
    }
}

void fill_matrix(fft_Matrix_cf* mat)
{
    for (size_t k = 0; k < mat->rows; ++k) {
        for (size_t j = 0; j < mat->cols; ++j) {
            if (k % 2 == 0) {
                fft_mat_set(mat, k, j, 13.f);
            } else {
                fft_mat_set(mat, k, j, 4.f);
            }
        }
    }
}

void save_img_as_ppm(uint32_t* img, size_t rows, size_t cols, const char *file_path)
{
    FILE* f = fopen(file_path, "wb");
    if (f == NULL) {
        fprintf(stderr, "ERROR: could not write into file %s", file_path);//: %s\n", file_path, strerror(errno));
        exit(1);
    }
    size_t img_w = cols;
    size_t img_h = rows;
    fprintf(f, "P6\n%ld %ld 255\n", img_w, img_h);
    for (size_t y=0; y<img_h; ++y) {
        for (size_t x=0; x<img_w; ++x) {
            uint32_t pixel = img[x+y*cols];
            uint8_t bytes[3] = {
                (pixel&0x0000FF)>>8*0,
                (pixel&0x00FF00)>>8*1,
                (pixel&0xFF0000)>>8*2
            };
            fwrite(bytes, sizeof(bytes), 1, f);
            assert(!ferror(f));
        }
    }
    fclose(f);
}

size_t sqr_dist(int x1, int y1, int x2, int y2)
{
    return (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
}

void draw_circle(int cx, int cy, int r, uint32_t* img, size_t img_w, size_t img_h) 
{
    for (size_t x = 0; x < img_w; ++x) {
        for (size_t y = 0; y < img_h; ++y) {
            if (sqr_dist(cx, cy, x, y) < r*r) {
                img[x+y*img_w] = 0xFF000000;
            } else {
                img[x+y*img_w] = 0xFFFFFFFF;
            }
        }
    }
}

void fill_image_white(uint32_t* img, size_t img_w, size_t img_h) 
{
    for (size_t x = 0; x < img_w; ++x) {
        for (size_t y = 0; y < img_h; ++y) {
            img[x+y*img_w] = 0xFFFFFFFF;
        }
    }
}

void rectangle(int rx, int ry, int w, int h, uint32_t* img, int img_w, int img_h) 
{
    for (int x = 0; x < img_w; ++x) {
        for (int y = 0; y < img_h; ++y) {
            if (x < rx+w && x > rx && y < ry+h && y > ry) {
                img[x+y*img_w] = 0xFF000000;
            }
        }
    }
}

void fill_image_sin_wave(uint32_t* img, int img_w, int img_h) 
{
    for (int x = 0; x < img_w; ++x) {
        for (int y = 0; y < img_h; ++y) {
            float coef = 200.f;
            float mag_x = (1.f + sin(x / coef)) / 2.f;
            float mag_y = (1.f + sin(y / coef)) / 2.f;
            uint32_t pixel = 0xFFFFFF * mag_x;
            pixel = (pixel + 0xFFFFFF * mag_y) / 2;
            img[x+y*img_w] = 0xFF000000 | pixel;
        }
    }
}

int main()
{
    // fft on the rows of a matrix
    // fft on the cols of the fft'd matrix
    size_t img_rows = 512, img_cols = 512;
    uint32_t* img = (uint32_t*)malloc(img_rows*img_cols*sizeof(uint32_t));
    // fill_image_white(img, img_cols, img_rows);
    // draw_circle(img_cols/2, img_rows/2, 15, img, img_cols, img_rows);
    // rectangle(img_cols/2 - 50, img_rows/2 - 50, 30, 10, img, img_cols, img_rows);
    fill_image_sin_wave(img, img_cols, img_rows);
    save_img_as_ppm(img, img_rows, img_cols, "original.ppm");

    fft_Matrix_cf mat = fft_mat_alloc(img_rows, img_cols);

    // size_t blr_rows = 4, blr_cols = 4;
    // fft_Matrix_cf blr = fft_mat_alloc(blr_rows, blr_cols);

    for (size_t x = 0; x < img_cols; ++x) {
        for (size_t y = 0; y < img_rows; ++y) {
            float norm = img[x+y*img_cols] / (float)0xFFFFFF;
            fft_mat_set(&mat, x, y, (float complex)norm);
        }
    }

    // for (size_t k = 0; k < blr.rows; ++k) {
        // for (size_t j = 0; j < blr.cols; ++j) {
            // fft_mat_set(&blr, k, j, 0.25f);
        // }
    // }

    printf("ORIGINAL:\n");
    printf("----------------------\n");
    // print_matrix(&mat);
    printf("----------------------\n");
    fft_matrix(&mat);
    // fft_matrix(&blr);
    for (size_t x = 0; x < img_cols; ++x) {
        for (size_t y = 0; y < img_rows; ++y) {
            float px = creal(fft_mat_get(&mat, x, y));
            // float sig = (1.f / 2.f*(1.f + exp(px))) + 1.f;
            uint32_t sig = px * 0xFFFFFF;
            // uint32_t mag = (sig * (float)0xFFFFFF);
            // printf("%f -> ", px);
            // printf("%f -> ", sig);
            // printf("%ld\n", mag);
            img[x+y*img_cols] = 0xFF000000|sig;
        }
    }
    save_img_as_ppm(img, img_rows, img_cols, "fft.ppm");

    // for (size_t k = 0; k < img_rows; ++k) {
        // for (size_t j = 0; j < img_cols; ++j) {
            // fft_mat_set(&mat, k, j, fft_mat_get(&mat, k, j) * fft_mat_get(&blr, k, j));
        // }
    // }
    // print_matrix(&mat);
    // printf("----------------------\n");
    // return 0;
    fft_matrix_inverse(&mat);
    for (size_t x = 0; x < img_cols; ++x) {
        for (size_t y = 0; y < img_rows; ++y) {
            float px = creal(fft_mat_get(&mat, x, y));
            // float sig = (1.f / (1.f + exp(px)));
            uint32_t mag = (px * 0xFFFFFF);
            // printf("%f -> ", px);
            // printf("%f -> ", sig);
            // printf("%ld\n", mag);
            img[x+y*img_cols] = 0xFF000000|mag;
        }
    }
    save_img_as_ppm(img, img_rows, img_cols, "restored.ppm");
    // printf("RESULT:\n");
    // printf("----------------------\n");
    // print_matrix(&mat);
    // printf("----------------------\n");

    fft_mat_free(&mat);
    // fft_mat_free(&blr);
    free(img);
    return 0;
}
