How does this compare to fc_consensus.py?

### Base speed
This is slightly faster on Centos7, but slower on Centos6.
Apparently, Centos6 has an odd inefficiency working with
threads, I guess.

### msa_array
The thread-local msa_array is lost when a thread is re-used.
(I do not know whether the thread is destroyed, but clearly thread-locals are lost.)
But memory is re-used from the heap pretty well now
(after araq fixed a bug in Nim).

### Memory
Memory is typically 550MB per thread, vs. 1GB per proc for Python.
So it is *faster* than Python for the same amount of memory
(2x as many procs).

### Missing funcs
We are missing with-trim functionality (though the func is
not called "without_trim"), as well as non-multiple output.
