## ##mangle uint64_t uint64

## #
## #  =====================================================================================
## # 
## #        Filename:  common.h
## # 
## #     Description:  Common delclaration for the code base 
## # 
## #         Version:  0.1
## #         Created:  07/16/2013 07:46:23 AM
## #        Revision:  none
## #        Compiler:  gcc
## # 
## #          Author:  Jason Chin, 
## #         Company:  
## # 
## #  =====================================================================================
## #
## # #################################################################################$$
## # # Copyright (c) 2011-2014, Pacific Biosciences of California, Inc.
## # #
## # # All rights reserved.
## # #
## # # Redistribution and use in source and binary forms, with or without
## # # modification, are permitted (subject to the limitations in the
## # # disclaimer below) provided that the following conditions are met:
## # #
## # #  * Redistributions of source code must retain the above copyright
## # #  notice, this list of conditions and the following disclaimer.
## # #
## # #  * Redistributions in binary form must reproduce the above
## # #  copyright notice, this list of conditions and the following
## # #  disclaimer in the documentation and/or other materials provided
## # #  with the distribution.
## # #
## # #  * Neither the name of Pacific Biosciences nor the names of its
## # #  contributors may be used to endorse or promote products derived
## # #  from this software without specific prior written permission.
## # #
## # # NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
## # # GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
## # # BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
## # # WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
## # # OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
## # # DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
## # # CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## # # SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## # # LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
## # # USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
## # # ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
## # # OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
## # # OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
## # # SUCH DAMAGE.
## # #################################################################################$$
## # 
from times import nil

type
  seq_coor_t* = cint
  alignment* = object
    aln_str_size*: seq_coor_t
    dist*: seq_coor_t
    aln_q_s*: seq_coor_t
    aln_q_e*: seq_coor_t
    aln_t_s*: seq_coor_t
    aln_t_e*: seq_coor_t
    q_aln_str*: string
    t_aln_str*: string

  d_path_data* = object
    pre_k*: seq_coor_t
    x1*: seq_coor_t
    y1*: seq_coor_t
    x2*: seq_coor_t
    y2*: seq_coor_t

  d_path_data2* = object
    d*: seq_coor_t
    k*: seq_coor_t
    pre_k*: seq_coor_t
    x1*: seq_coor_t
    y1*: seq_coor_t
    x2*: seq_coor_t
    y2*: seq_coor_t

  path_point* = object
    x*: seq_coor_t
    y*: seq_coor_t

  kmer_lookup* = object
    start*: seq_coor_t
    last*: seq_coor_t
    count*: seq_coor_t

  base* = cuchar
  seq_array* = seq[base] #ptr base
  seq_addr* = seq_coor_t
  seq_addr_array* = ptr seq_addr
  kmer_match* = object
    count*: seq_coor_t # triply redunant
    query_pos*: seq[seq_coor_t]
    target_pos*: seq[seq_coor_t]

  aln_range* = object
    s1*: seq_coor_t
    e1*: seq_coor_t
    s2*: seq_coor_t
    e2*: seq_coor_t
    score*: clong

  consensus_data* = object
    sequence*: string
    eqv*: seq[cint]

# For ptr arithmetic
template usePtr*[T] =
  template `+`(p: ptr T, off: Natural): ptr T =
    cast[ptr type(p[])](cast[ByteAddress](p) +% int(off) * sizeof(p[]))

  template `+=`(p: ptr T, off: Natural) =
    p = p + off

  template `-`(p: ptr T, off: Natural): ptr T =
    cast[ptr type(p[])](cast[ByteAddress](p) -% int(off) * sizeof(p[]))

  template `-=`(p: ptr T, off: Natural) =
    p = p - int(off)

  template `[]`(p: ptr T, off: Natural): T =
    (p + int(off))[]

  template `[]=`(p: ptr T, off: Natural, val: T) =
    (p + off)[] = val

# https://forum.nim-lang.org/t/1188/1
template ptrMath*(body: untyped) =
  template `+`[T](p: ptr T, off: Natural): ptr T =
    cast[ptr type(p[])](cast[ByteAddress](p) +% int(off) * sizeof(p[]))
  
  template `+=`[T](p: ptr T, off: Natural) =
    p = p + off
  
  template `-`[T](p: ptr T, off: Natural): ptr T =
    cast[ptr type(p[])](cast[ByteAddress](p) -% int(off) * sizeof(p[]))
  
  template `-=`[T](p: ptr T, off: Natural) =
    p = p - int(off)
  
  template `[]`[T](p: ptr T, off: Natural): T =
    (p + int(off))[]
  
  template `[]=`[T](p: ptr T, off: Natural, val: T) =
    (p + off)[] = val
  
  body

proc calloc*[T](n: Natural): ptr T =
  return cast[ptr T](alloc(n * sizeof(T)))
proc calloc0*[T](n: Natural): ptr T =
  return cast[ptr T](alloc0(n * sizeof(T)))
proc realloc*[T](p: ptr T, n: Natural): ptr T =
  return cast[ptr T](alloc(n * sizeof(T)))
proc realloc0*[T](p: ptr T, n: Natural): ptr T =
  return cast[ptr T](alloc0(n * sizeof(T)))

proc log*(msgs: varargs[string]) =
  for s in msgs:
    write(stderr, s)
  write(stderr, '\l')
  return
template benchmark*(benchmarkName: string, code: typed) =
  let t0 = times.epochTime()
  code
  let elapsed = times.epochTime() - t0
  let elapsedStr = elapsed.formatFloat(format = ffDecimal, precision = 3)
  log("CPU Time [$1] $2s" % [$benchmarkName, $elapsedStr])
