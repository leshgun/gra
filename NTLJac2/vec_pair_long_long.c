#include "vec_pair_long_long.h"

#include <NTL/new.h>

NTL_START_IMPL

NTL_pair_impl(long,long,pair_long_long)
NTL_pair_io_impl(long,long,pair_long_long)
NTL_pair_eq_impl(long,long,pair_long_long)

NTL_vector_impl(pair_long_long,vec_pair_long_long)
NTL_io_vector_impl(pair_long_long,vec_pair_long_long)
NTL_eq_vector_impl(pair_long_long,vec_pair_long_long)

NTL_END_IMPL

