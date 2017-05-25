#include "data_sorter.h"
#include <assert.h>
#include <stdlib.h>

int compare_path_data2(const void* a, const void* b)
{
    const d_path_data2* arg1 = a;
    const d_path_data2* arg2 = b;
    if (arg1->d - arg2->d == 0) {
        return  arg1->k - arg2->k;
    } else {
        return arg1->d - arg2->d;
    }
}


void d_path_data2sort(d_path_data2* base, size_t max_idx, size_t size) {
    assert(size == sizeof(d_path_data2));
    qsort(base, max_idx, size, compare_path_data2);
}

void do_new_x_and_y(seq_coor_t* x0, seq_coor_t* y0, seq_coor_t q_len, seq_coor_t t_len,
        char const* query_seq, char const* target_seq)
{
    seq_coor_t x = *x0;
    seq_coor_t y = *y0;
            while ( x < q_len && y < t_len && query_seq[x] == target_seq[y] ){
                x++;
                y++;
            }
    *x0 = x; *y0 = y;
}
