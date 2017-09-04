#include <stddef.h>

typedef int seq_coor_t; 

typedef struct {
    seq_coor_t d;
    seq_coor_t k;
    seq_coor_t pre_k;
    seq_coor_t x1;
    seq_coor_t y1;
    seq_coor_t x2;
    seq_coor_t y2;
} d_path_data2;

extern void d_path_data2sort(d_path_data2* base, size_t length, size_t max_idx);

extern void do_new_x_and_y(seq_coor_t* x0, seq_coor_t* y0, seq_coor_t q_len, seq_coor_t t_len,
        char const* query_seq, char const* target_seq);
