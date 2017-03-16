The problem with calling the Falcon C code is that the msa_array is global.

If we turn off STATIC_ALLOCATE, then things are very slow.

Anyway, this might be interesting for future reference.
