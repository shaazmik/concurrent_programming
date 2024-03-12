#include <math.h>
#include <stdint.h>
#include <stdio.h>

typedef struct dim {
    double len;
    uint32_t pts_num;
    double step;
} dim;

typedef struct dimensions {
    dim x;
    dim t;
} dimensions;

typedef struct data_xy {
    double **data;
    dimensions d;
    double d_coeff;
} data_xy;



static double **data_alloc(uint32_t pts_x, uint32_t pts_y) {
    double **data = (double **) calloc(pts_x, sizeof(data[0]));
    
    for (uint32_t i = 0; i < pts_x; i++) {
        data[i] = (double *) calloc(pts_y, sizeof(data[i][0]));
    }

    return data;
}

static void data_print(data_xy *state) {
    for (uint32_t i = 0; i < state->d.x.pts_num; i++) {
        printf("%g ", state->data[i][0]);
        putc('\n', stdout);
    }
}

static double calculate_new_point(data_xy *state, uint32_t x_pos, uint32_t y_pos) {
    double **data = state->data;

    double r_side = (data[x_pos + 1][y_pos] - 2 * data[x_pos][y_pos] + data[x_pos - 1][y_pos]) / state->d.x.step / state->d.x.step;

    return data[x_pos][y_pos] + state->d.t.step * state->d_coeff * r_side;
}

static void calculate_row(data_xy *state, uint32_t row) {
    if ((row < 1) || (row >= state->d.t.pts_num)) return;
    for (uint32_t i = 1; i < state->d.x.pts_num - 1; i++) {
        state->data[i][row] = calculate_new_point(state, i, row - 1);
    }
}

static void data_free(data_xy *state) {
    for (uint32_t i = 0; i < state->d.x.pts_num; i++) {
        free(state->data[i]);
    }

    free(state->data);
}

static void init_start_coord_dim(double** data, uint32_t pts_num, double step) {
    for (uint32_t i = 0; i < pts_num; i++) {
        data[i][0] = sin(step * i * 5);
    }
}

static void init_start_time_dim(double** data, uint32_t pts_num, double step) {
    for (uint32_t i = 0; i < pts_num; i++) {
        data[0][i] = data[0][0] + (step * i);
    }
}


int main(int argc, const char **argv) {
    dim x = {1.0, 1000, 0.0};
    x.step = x.len / x.pts_num;

    dim t = {20.0, 1000000, 0.0};
    t.step = t.len / t.pts_num;

    const dimensions dims = {x, t};

    double **data = data_alloc(x.pts_num, t.pts_num);

    init_start_coord_dim(data, x.pts_num, x.step);

    init_start_time_dim(data, t.pts_num, t.step);

    // initial conditions at the end of the grid for the time space
    for (uint32_t i = 1; i < t.pts_num; i++) {
        data[x.pts_num - 1][i] = data[x.pts_num - 1][0] - (t.step * i);
    }

    data_xy state = {data, dims, 1e-3};

    for (uint32_t row = 1; row < t.pts_num; row++) {
        calculate_row(&state, row);
    }

    for (uint32_t i = 1; i < x.pts_num; i++) {
        printf("data[%d][%d] = %lf\n", i, t.pts_num - 1, data[i][t.pts_num - 1]);
    }
    data_free(&state);

    return 0;
}
