## vim: sw=2 ts=2 sts=2 tw=0 et:
## #
## #  =====================================================================================
## # 
## #        Filename:  DW_banded.c
## # 
## #     Description:  A banded version for the O(ND) greedy sequence alignment algorithm
## # 
## #         Version:  0.1
## #         Created:  07/20/2013 17:00:00
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
## #
## #

{.passC: "-Wall -I.".}
#{.passC: "-g -Wall -I. -fno-omit-frame-pointer -fno-strict-aliasing -DNDEBUG".} // from pip install -v
{.compile: "data_sorter.c".}

import common
#import algorithm  # We need a sort with a max-index.

common.usePtr[byte]()
common.usePtr[char]()
#common.usePtr[d_path_data2]()

# Was used by sort. Still used by bsearch.
proc compare_d_path*(arg1, arg2: d_path_data2): int =
  #log("Comparing:", repr(arg1), " to ", repr(arg2))
  if arg1.d == arg2.d:
    if arg1.k == arg2.k:
      result = 0
    else:
      result = if arg1.k < arg2.k: -1 else: 1
  else:
    result = if arg1.d < arg2.d: -1 else: 1
  #log("result=", $result)

proc d_path_data2sort*(data: pointer; max_idx: ByteAddress, size: ByteAddress) {.
    cdecl, importc: "d_path_data2sort", header: "data_sorter.h".}

proc d_path_sort*(path_base: var seq[d_path_data2], max_idx: int32) =
  if max_idx == 1:
    return
  #log("sort:", $len(path_base), " max_idx:", $max_idx)
  #path_base.sort(compare_d_path)
  d_path_data2sort(pointer(addr path_base[0]), ByteAddress(max_idx), ByteAddress(sizeof(d_path_data2)));
  #void d_path_sort( d_path_data2 * base, unsigned long max_idx)
  #qsort(base, max_idx, sizeof(d_path_data2), compare_d_path);

# Copied from Nim/lib/pure/algorithm.nim, and added cmp proc,
# and max_idx.
proc binarySearch*[T](a: openArray[T], length: int, key: T, cmp: proc (x, y: T): int {.closure.}): int =
  #log("length:", $length, ", key:", $key)
  ## binary search for `key` in `a`. Returns -1 if not found.
  var b = length # must be <= len(a)
  result = 0
  while result < b:
    var mid = (result + b) div 2
    if cmp(a[mid], key) < 0: result = mid + 1
    else: b = mid
  if result >= length or cmp(a[result], key) != 0: result = -1

proc get_dpath_idx(d: seq_coor_t; k: seq_coor_t; max_idx: int; # culong?
                   path_base: seq[d_path_data2]): d_path_data2 =
  var rtn: d_path_data2
  var d_tmp: d_path_data2
  d_tmp.d = d
  d_tmp.k = k
  #log("d_tmp:", repr(d_tmp), ", len(path_base)=", $len(path_base), ", max_idx=", $max_idx)
  #log("path_base[0]:", repr(path_base[0]))
  #log("path_base[1]:", repr(path_base[1]))
  #log(repr(path_base[2236]))
  #log("...")
  var found = binarySearch(path_base, max_idx, d_tmp, compare_d_path)
  #log("found:", $found)
  if found == -1:
    raise newException(ValueError, "Could not find d,k in path_base") # TODO(CD): encode key into msg
  rtn = path_base[found]
  #log("dp ", $repr(rtn))
  ## #printf("dp %ld %ld %ld %ld %ld %ld %ld\n", (rtn)->d, (rtn)->k, (rtn)->x1, (rtn)->y1, (rtn)->x2, (rtn)->y2, (rtn)->pre_k);
  return rtn

#discard """
proc print_d_path*(base: var seq[d_path_data2]; max_idx: int) =
  var idx: int
  idx = 0
  while idx < max_idx:
    echo "dp ", idx, " ", repr(base[idx.int])
    inc(idx)
#"""

proc bzero[T](uninit: var seq[T]) =
  let n = len(uninit) * sizeof(T)
  let bz: ptr byte = cast[ptr byte](addr (uninit[0]))
  system.zeroMem(bz, n)
proc bzero(uninit: var string) =
  let n = len(uninit)
  let bz: ptr char = cast[ptr char](addr (uninit[0]))
  zeroMem(bz, n)

#
# These are split from align() only for profiling.
#

proc mem(d_path: var seq[d_path_data2], max_d, band_size: seq_coor_t) =
  ## # We should probably use hashmap to store the backtracing information to save memory allocation time
  ## # This O(MN) block allocation scheme is convient for now but it is slower for very long sequences
  let ssize = max_d * (band_size + 1) * 2 + 1
  #if ssize > 1000000:
  #  raise newException(ValueError, "too big") # Just to catch bugs during development.

  #log("Big seq:", $ssize)
  #newSeq(d_path, (max_d * (band_size + 1) * 2 + 1)) # maybe drop +1?
  d_path.setLen(ssize)
  #bzero(d_path)
  ## #fprintf(stderr, "calloc(%d x %d)\n", max_d * (band_size + 1 ) * 2 + 1, sizeof(d_path_data2));
proc mem2(aln_path: var seq[path_point], q_len, t_len: seq_coor_t, align_rtn: var alignment) =
  #var align_rtn: ref alignment
  #newSeq(aln_path, (q_len + t_len + 1))
  #new(align_rtn)
  aln_path.setLen(q_len + t_len + 1)
  #bzero(aln_path)
  let slen = q_len + t_len + 1
  align_rtn.t_aln_str.setLen(slen) # = newString(q_len + t_len + 1)
  align_rtn.q_aln_str.setLen(slen) # = newString(q_len + t_len + 1)
  #bzero(align_rtn.t_aln_str)
  #bzero(align_rtn.q_aln_str)
  align_rtn.aln_str_size = 0
  align_rtn.aln_q_s = 0
  align_rtn.aln_q_e = 0
  align_rtn.aln_t_s = 0
  align_rtn.aln_t_e = 0
  # return align_rtn
proc do_get_aln_str(d, k, q_len, t_len: seq_coor_t, max_idx: int32, target_seq, query_seq: ptr char, d_path: var seq[d_path_data2], aln_path: var seq[path_point], align_rtn: var alignment) =
        var cd: seq_coor_t
        var ck: seq_coor_t
        var aln_path_idx: seq_coor_t
        var
          cx: seq_coor_t
          cy: seq_coor_t
          nx: seq_coor_t
          ny: seq_coor_t
        cd = d
        ck = k
        aln_path_idx = 0
        while cd >= 0 and aln_path_idx < q_len + t_len + 1:
          let d_path_aux: d_path_data2 = get_dpath_idx(cd, ck, max_idx, d_path)
          #log("d_path_aux:",repr(d_path_aux))
          #log("aln_path_idx:", $aln_path_idx)
          #log("len:", $len(aln_path))
          #log("sub:", repr(aln_path[aln_path_idx]))
          aln_path[aln_path_idx].x = d_path_aux.x2
          aln_path[aln_path_idx].y = d_path_aux.y2
          inc(aln_path_idx)
          aln_path[aln_path_idx].x = d_path_aux.x1
          aln_path[aln_path_idx].y = d_path_aux.y1
          inc(aln_path_idx)
          ck = d_path_aux.pre_k
          dec(cd, 1)
        dec(aln_path_idx)
        cx = aln_path[aln_path_idx].x
        cy = aln_path[aln_path_idx].y
        align_rtn.aln_q_s = cx
        align_rtn.aln_t_s = cy
        var aln_pos: seq_coor_t
        aln_pos = 0
        while aln_path_idx > 0:
          dec(aln_path_idx)
          nx = aln_path[aln_path_idx].x
          ny = aln_path[aln_path_idx].y
          if cx == nx and cy == ny:
            continue
          var i: seq_coor_t
          if nx == cx and ny != cy:
            ## #advance in y
            i = 0
            while i < ny - cy:
              align_rtn.q_aln_str[aln_pos + i] = '-'
              inc(i)
            i = 0
            while i < ny - cy:
              align_rtn.t_aln_str[aln_pos + i] = target_seq[cy + i]
              inc(i)
            inc(aln_pos, ny - cy)
          elif nx != cx and ny == cy:
            ## #advance in x
            i = 0
            while i < nx - cx:
              align_rtn.q_aln_str[aln_pos + i] = query_seq[cx + i]
              inc(i)
            i = 0
            while i < nx - cx:
              align_rtn.t_aln_str[aln_pos + i] = '-'
              inc(i)
            inc(aln_pos, nx - cx)
          else:
            i = 0
            while i < nx - cx:
              align_rtn.q_aln_str[aln_pos + i] = query_seq[cx + i]
              inc(i)
            i = 0
            while i < ny - cy:
              align_rtn.t_aln_str[aln_pos + i] = target_seq[cy + i]
              inc(i)
            inc(aln_pos, ny - cy)
          cx = nx
          cy = ny
        align_rtn.aln_str_size = aln_pos
proc new_min_and_max(U: seq[seq_coor_t], min_k, max_k, k_offset, best_m, band_tolerance: seq_coor_t): (seq_coor_t, seq_coor_t) =
    ## # For banding
    var
      new_min_k = max_k
      new_max_k = min_k
    var k2 = min_k
    while k2 <= max_k:
      if U[k2 + k_offset] >= best_m - band_tolerance:
        if k2 < new_min_k:
          new_min_k = k2
        if k2 > new_max_k:
          new_max_k = k2
      inc(k2, 2)
    return (new_min_k - 1, new_max_k + 1)
# Alternative C version.
proc do_new_x_and_y*(x0, y0: ptr seq_coor_t, q_len, t_len: seq_coor_t, query_seq, target_seq: ptr char) {.
    cdecl, importc: "do_new_x_and_y", header: "data_sorter.h".}
proc new_x_and_y(x0, y0, q_len, t_len: seq_coor_t, query_seq, target_seq: ptr char): auto =
      var x = x0
      var y = y0
      while x < q_len and y < t_len and query_seq[x] == target_seq[y]:
        inc(x)
        inc(y)
      return (x, y)
proc align1*(query_seq: ptr char; q_len: seq_coor_t; target_seq: ptr char;
           t_len: seq_coor_t; band_tolerance: seq_coor_t; get_aln_str: bool,
           d_path: var seq[d_path_data2], aln_path: var seq[path_point], align_rtn: var alignment, U,V: var seq[seq_coor_t]) =
  #log("In align...")
  ## # array of matched bases for each "k"
  var k_offset: seq_coor_t
  var d: seq_coor_t
  var
    k: seq_coor_t
  var best_m: seq_coor_t
  ## # the best "matches" for each d
  var
    min_k: seq_coor_t
  var
    max_k: seq_coor_t
  var pre_k: seq_coor_t
  var
    x: seq_coor_t
    y: seq_coor_t
  var max_d: seq_coor_t
  var band_size: seq_coor_t
  var d_path_idx: int32 = 0
  var max_idx: int32 = 0
  #var d_path: seq[d_path_data2]
  #var aln_path: seq[path_point]
  #var align_rtn: ref alignment
  var aligned: bool = false
  ## #printf("debug: %ld %ld\n", q_len, t_len);
  ## #printf("%s\n", query_seq);
  max_d = seq_coor_t(0.3 * float(q_len + t_len))
  let uv_len = max_d * 2 + 1
  U.setLen(uv_len)
  V.setLen(uv_len)
  # These bzero() are needed, unfortunately. Can they be faster?
  bzero(U)
  bzero(V)
  band_size = band_tolerance * 2
  mem(d_path, max_d, band_size)
  mem2(aln_path, q_len, t_len, align_rtn)
  ## #fprintf(stderr, "UV sz aln sz2: %d %d %d %d\n", (max_d * 2 + 1), sizeof(seq_coor_t), (q_len + t_len + 1), sizeof(d_path_data2));
  ## #printf("max_d: %lu, band_size: %lu\n", max_d, band_size);
  k_offset = max_d
  best_m = - 1
  min_k = 0
  max_k = 0
  d_path_idx = 0
  max_idx = 0
  d = 0
  while d < max_d:
    if max_k - min_k > band_size:
      break
    k = min_k
    while k <= max_k:
      if (k == min_k) or
          ((k != max_k) and (V[k - 1 + k_offset] < V[k + 1 + k_offset])):
        pre_k = k + 1
        x = V[k + 1 + k_offset]
      else:
        pre_k = k - 1
        x = V[k - 1 + k_offset] + 1
      y = x - k
      d_path[d_path_idx].d = d
      d_path[d_path_idx].k = k
      d_path[d_path_idx].x1 = x
      d_path[d_path_idx].y1 = y
      (x, y) = new_x_and_y(x, y, q_len, t_len, query_seq, target_seq)
      # This takes the most time, according to profiler. I am skeptical, but
      # maybe the C version can be faster somehow.
      #do_new_x_and_y(addr x, addr y, q_len, t_len, query_seq, target_seq)
      d_path[d_path_idx].x2 = x
      d_path[d_path_idx].y2 = y
      d_path[d_path_idx].pre_k = pre_k
      inc(d_path_idx)
      V[k + k_offset] = x
      U[k + k_offset] = x + y
      if x + y > best_m:
        best_m = x + y
      if x >= q_len or y >= t_len:
        aligned = true
        max_idx = d_path_idx
        break
      inc(k, 2)
    ## # For banding
    (min_k, max_k) = new_min_and_max(U, min_k, max_k, k_offset, best_m, band_tolerance)
    ## # For no banding
    ## # max_k ++;
    ## # min_k --;
    ## # For debuging
    #printf("min_max_k,d, %ld %ld %ld\n", min_k, max_k, d);
    if aligned == true:
      align_rtn.aln_q_e = x
      align_rtn.aln_t_e = y
      align_rtn.dist = d
      align_rtn.aln_str_size = (x + y + d) div 2
      align_rtn.aln_q_s = 0
      align_rtn.aln_t_s = 0
      #print_d_path(d_path, max_idx);
      #echo "================="
      d_path_sort(d_path, max_idx)
      #print_d_path(d_path, max_idx);
      #log("align_rtn0:", $align_rtn[])
      if get_aln_str:
        do_get_aln_str(d, k, q_len, t_len, max_idx, target_seq, query_seq, d_path, aln_path, align_rtn)
      break
    inc(d)
  #free(aln_path)
  #log("Leaving align")
  #return align_rtn

proc align*(query_seq: ptr char; q_len: seq_coor_t; target_seq: ptr char;
           t_len: seq_coor_t; band_tolerance: seq_coor_t; get_aln_str: bool): ref alignment =
  var align_rtn: ref alignment
  new(align_rtn)
  var d_path: seq[d_path_data2]
  var aln_path: seq[path_point]
  var U: seq[seq_coor_t]
  var V: seq[seq_coor_t]
  newSeq(d_path, 0)
  newSeq(aln_path, 0)
  align_rtn.t_aln_str = newString(0)
  align_rtn.q_aln_str = newString(0)
  newSeq(U, 0)
  newSeq(V, 0)
  align1(query_seq, q_len, target_seq, t_len, band_tolerance, get_aln_str,
      d_path, aln_path, align_rtn[], U, V)
  return align_rtn

#proc free_alignment*(aln: ptr alignment) =
#  free(aln.q_aln_str)
#  free(aln.t_aln_str)
#  free(aln)
