# vim: sw=2 ts=2 sts=2 tw=80 et:
## #
## #  =====================================================================================
## # 
## #        Filename:  kmer_count.c
## # 
## #     Description:
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

import common
import algorithm

common.usePtr[kmer_lookup]
common.usePtr[common.base]
common.usePtr[common.seq_addr]

const
  UINT8_MAX = 255
  UINT16_MAX = uint16.high.int32
  UINT_MAX = uint32.high
  INT_MAX = int32.high
  LONG_MIN = clong.low

var KMERMATCHINC*: seq_coor_t = 10000
discard """
proc compare_seq_coor*(a: pointer; b: pointer): cint =
  var arg1: ptr seq_coor_t = a
  var arg2: ptr seq_coor_t = b
  return (arg1[]) - (arg2[])
"""
proc init_kmer_lookup*(kl: var seq[kmer_lookup]; size: seq_coor_t) =
  var i: seq_coor_t
  ## #printf("%lu is allocated for kmer lookup\n", size);
  i = 0
  while i < size:
    kl[i].start = INT_MAX
    kl[i].last = INT_MAX
    kl[i].count = 0
    inc(i)

proc allocate_kmer_lookup*(size: seq_coor_t): seq[kmer_lookup] =
  newSeq(result, size)
  ## #printf("%lu is allocated for kmer lookup\n", size);
  init_kmer_lookup(result, size)

proc init_seq_array*(sa: var seq_array; size: seq_coor_t) =
  var i: seq_coor_t
  i = 0
  while i < size:
    sa[i] = 0x000000FF.base
    inc(i)

proc allocate_seq*(size: seq_coor_t): seq_array =
  newSeq(result, size)
  init_seq_array(result, size)

proc allocate_seq_addr*(size: seq_coor_t): seq_addr_array =
  return calloc0[seq_addr](size)
  # Apparently, find_kmer_pos_for_seq() loop requires 0s

proc get_kmer_bitvector*(sa: ptr base; K: cuint): seq_coor_t =
  var i: cuint
  var kmer_bv: seq_coor_t = 0
  var kmer_mask: seq_coor_t
  kmer_mask = 0 # TODO(CD): Why is kmer_mask completely ignored??
  i = 0
  while i < K:
    kmer_mask = kmer_mask shl 2
    kmer_mask = kmer_mask or 0x00000003
    inc(i)
  i = 0
  while i < K:
    kmer_bv = kmer_bv shl 2
    kmer_bv = kmer_bv or seq_coor_t((int(sa[i.int])) and 0x00000003)
    inc(i)
  return kmer_bv

proc add_sequence*(start: seq_coor_t; K: cuint; cseq: cstring; seq_len: seq_coor_t;
                  sda: seq_addr_array; sa: var seq_array; lk: var seq[kmer_lookup]) =
  var i: seq_coor_t
  var kmer_bv: seq_coor_t
  var kmer_mask: seq_coor_t
  kmer_mask = 0
  i = 0
  #log("Adding seq of length:", $seq_len)
  while i < K.seq_coor_t:
    kmer_mask = kmer_mask shl 2
    kmer_mask = kmer_mask or 0x00000003
    inc(i)
  i = 0
  while i < seq_len:
    case cseq[i]
    of 'A':
      sa[start + i] = 0.base
    of 'C':
      sa[start + i] = 1.base
    of 'G':
      sa[start + i] = 2.base
    of 'T':
      sa[start + i] = 3.base
    else:
       raise newException(ValueError, "Must be ACGT")
    inc(i)
  kmer_bv = get_kmer_bitvector(addr sa[start], K)
  i = 0
  while i < seq_len - K.seq_coor_t:
    ## #printf("%lu %lu\n", i, kmer_bv);
    ## #printf("lk before init: %lu %lu %lu\n", kmer_bv, lk[kmer_bv].start, lk[kmer_bv].last);
    if lk[kmer_bv].start == INT_MAX:
      lk[kmer_bv].start = start + i
      lk[kmer_bv].last = start + i
      inc(lk[kmer_bv].count, 1)
      ## #printf("lk init: %lu %lu %lu\n", kmer_bv, lk[kmer_bv].start, lk[kmer_bv].last);
    else:
      sda[lk[kmer_bv].last] = seq_addr(start + i)
      inc(lk[kmer_bv].count, 1)
      lk[kmer_bv].last = start + i
      ## #printf("lk change: %lu %lu %lu\n", kmer_bv, lk[kmer_bv].start, lk[kmer_bv].last);
    kmer_bv = kmer_bv shl 2
    kmer_bv = kmer_bv or seq_coor_t(sa[start + i + K.int32])
    kmer_bv = kmer_bv and kmer_mask
    inc(i)

proc mask_k_mer*(size: seq_coor_t; kl: ptr kmer_lookup; threshold: seq_coor_t) =
  var i: seq_coor_t
  i = 0
  while i < size:
    if kl[i].count > threshold:
      kl[i].start = INT_MAX
      kl[i].last = INT_MAX
      ## #kl[i].count = 0;
    inc(i)

proc find_kmer_pos_for_seq*(cseq: cstring; seq_len: seq_coor_t; K: cuint;
                           sda: seq_addr_array; lk: seq[kmer_lookup]): ref kmer_match =
  var i: seq_coor_t
  var kmer_bv: seq_coor_t
  var kmer_mask: seq_coor_t
  var kmer_pos: seq_coor_t
  var next_kmer_pos: seq_coor_t
  var half_K: cuint
  #var result_allocation_size: seq_coor_t = KMERMATCHINC
  var sa: seq[base]
  new(result) # was kmer_match_rtn
  result.count = 0 # redundant now, since the seqs know their lengths
  newSeq(result.query_pos, 0)
  newSeq(result.target_pos, 0)
  newSeq(sa, seq_len)
  kmer_mask = 0
  i = 0
  while i < K.seq_coor_t:
    kmer_mask = kmer_mask shl 2
    kmer_mask = kmer_mask or 0x00000003
    inc(i)
  i = 0
  while i < seq_len:
    case cseq[i]
    of 'A':
      sa[i] = 0.base
    of 'C':
      sa[i] = 1.base
    of 'G':
      sa[i] = 2.base
    of 'T':
      sa[i] = 3.base
    else:
      raise newException(ValueError, "Must be ACGT for kmer")
    inc(i)
  # TODO(CD): Isn't this call redundant later?
  kmer_bv = get_kmer_bitvector(addr sa[0], K)
  half_K = K shr 1
  #echo "hk:", half_K, " K:", K, " bv:", repr(kmer_bv)
  i = 0
  while i < seq_len - K.seq_coor_t:
    kmer_bv = get_kmer_bitvector(addr sa[i], K)
    if lk[kmer_bv].start == INT_MAX:
      ## #for high count k-mers
      inc(i, half_K.seq_coor_t)
      continue
    kmer_pos = lk[kmer_bv].start
    next_kmer_pos = sda[kmer_pos]
    result.query_pos.add(i)
    result.target_pos.add(kmer_pos)
    inc(result.count, 1)
    while next_kmer_pos > kmer_pos:
      kmer_pos = next_kmer_pos
      #echo "MY SDA:", repr(sda), " kmer_pos:", kmer_pos, " =", sda[kmer_pos]
      next_kmer_pos = sda[kmer_pos]
      result.query_pos.add(i)
      result.target_pos.add(kmer_pos)
      inc(result.count, 1)
    inc(i, half_K.seq_coor_t)
  #log("kmer_match:@", $kmer_pos, " ", $result.count, " ", $len(result.query_pos))
  #return result # implicit

proc find_best_aln_range*(km_ptr: ref kmer_match; bin_size: int;
                         count_th: seq_coor_t): ref aln_range =
  #log("find_best_aln_range(", $bin_size, ", ", $count_th, ")")
  var i: seq_coor_t
  var j: seq_coor_t
  var
    q_min: seq_coor_t
    q_max: seq_coor_t
    t_min: seq_coor_t
    t_max: seq_coor_t
  var d_count: seq[seq_coor_t]
  var q_coor: seq[seq_coor_t]
  var t_coor: seq[seq_coor_t]
  var arange: ref aln_range
  var
    d: clong
    d_min: clong
    d_max: clong
  var cur_score: clong
  var max_score: clong
  var max_k_mer_count: clong
  var max_k_mer_bin: clong
  var cur_start: seq_coor_t
  new(arange)
  q_min = INT_MAX
  q_max = 0
  t_min = INT_MAX
  t_max = 0
  d_min = INT_MAX
  d_max = LONG_MIN
  i = 0
  while i < km_ptr.count:
    if km_ptr.query_pos[i] < q_min:
      q_min = km_ptr.query_pos[i]
    if km_ptr.query_pos[i] > q_max:
      q_max = km_ptr.query_pos[i]
    if km_ptr.target_pos[i] < t_min:
      t_min = km_ptr.target_pos[i]
    if km_ptr.query_pos[i] > t_max:
      t_max = km_ptr.target_pos[i]
    d = clong(km_ptr.query_pos[i]) - clong(km_ptr.target_pos[i])
    if d < d_min:
      d_min = d
    if d > d_max:
      d_max = d
    inc(i)
  #log("KMM:", $km_ptr.count, " ", $d_min, " ", $d_max, " ", $((d_max - d_min) div bin_size + 1))
  newSeq(d_count, ((d_max - d_min) div bin_size + 1))
  newSeq(q_coor, km_ptr.count)
  newSeq(t_coor, km_ptr.count)
  i = 0
  while i < km_ptr.count:
    d = clong((km_ptr.query_pos[i])) - clong((km_ptr.target_pos[i]))
    inc(d_count[(d - d_min) div clong(bin_size)], 1)
    q_coor[i] = INT_MAX
    t_coor[i] = INT_MAX
    inc(i)
  j = 0
  max_k_mer_count = 0
  max_k_mer_bin = INT_MAX
  i = 0
  while i < km_ptr.count:
    d = clong((km_ptr.query_pos[i])) - clong((km_ptr.target_pos[i]))
    if d_count[(d - d_min) div clong(bin_size)] > max_k_mer_count:
      max_k_mer_count = d_count[(d - d_min) div clong(bin_size)]
      max_k_mer_bin = (d - d_min) div clong(bin_size)
    inc(i)
  #log("k_mer:", $max_k_mer_count, " ", $max_k_mer_bin)
  if max_k_mer_bin != INT_MAX and max_k_mer_count > count_th:
    i = 0
    while i < km_ptr.count:
      d = clong((km_ptr.query_pos[i])) - clong((km_ptr.target_pos[i]))
      if abs(((d - d_min) div clong(bin_size)) - max_k_mer_bin) > 5:
        inc(i)
        continue
      if d_count[(d - d_min) div clong(bin_size)] > count_th:
        q_coor[j] = km_ptr.query_pos[i]
        t_coor[j] = km_ptr.target_pos[i]
        ## #printf("d_count: %lu %lu\n" ,i, d_count[(d - d_min)/ (long int) bin_size]);
        ## #printf("coor: %lu %lu\n" , q_coor[j], t_coor[j]);
        inc(j)
      inc(i)
  if j > 1:
    arange.s1 = q_coor[0]
    arange.e1 = q_coor[0]
    arange.s2 = t_coor[0]
    arange.e2 = t_coor[0]
    arange.score = 0
    max_score = 0
    cur_score = 0
    cur_start = 0
    i = 1
    while i < j:
      inc(cur_score, 32 - (q_coor[i] - q_coor[i - 1]))
      #log("deltaD:", $(q_coor[i] - q_coor[i-1]), " ", $cur_score)
      if cur_score < 0:
        cur_score = 0
        cur_start = i
      elif cur_score > max_score:
        arange.s1 = q_coor[cur_start]
        arange.s2 = t_coor[cur_start]
        arange.e1 = q_coor[i]
        arange.e2 = t_coor[i]
        max_score = cur_score
        arange.score = max_score
        #log("Arange:", repr(arange))
      inc(i)
  else:
    arange.s1 = 0
    arange.e1 = 0
    arange.s2 = 0
    arange.e2 = 0
    arange.score = 0
  ## # printf("free\n");
  return arange

proc find_best_aln_range2*(km_ptr: ptr kmer_match;
                          bin_width: seq_coor_t; count_th: seq_coor_t): ref aln_range =
  var d_coor: seq[seq_coor_t]
  var hit_score: seq[seq_coor_t]
  var hit_count: seq[seq_coor_t]
  var last_hit: seq[seq_coor_t]
  var
    max_q: seq_coor_t
    max_t: seq_coor_t
  var
    s: seq_coor_t
    e: seq_coor_t
    max_s: seq_coor_t
    max_e: seq_coor_t
    max_span: seq_coor_t
    d_s: seq_coor_t
    d_e: seq_coor_t
    delta: seq_coor_t
    d_len: seq_coor_t
  var
    px: seq_coor_t
    py: seq_coor_t
    cx: seq_coor_t
    cy: seq_coor_t
  var max_hit_idx: seq_coor_t
  var
    max_hit_score: seq_coor_t
    max_hit_count: seq_coor_t = 0
  var
    i: seq_coor_t
    j: seq_coor_t
  var
    candidate_idx: seq_coor_t
    max_d: seq_coor_t
    d: seq_coor_t
  var arange: ref aln_range
  new(arange)
  newSeq(d_coor, km_ptr.count)
  max_q = - 1
  max_t = - 1
  i = 0
  while i < km_ptr.count:
    d_coor[i] = km_ptr.query_pos[i] - km_ptr.target_pos[i]
    max_q = if max_q > km_ptr.query_pos[i]: max_q else: km_ptr.query_pos[i]
    max_t = if max_t > km_ptr.target_pos[i]: max_q else: km_ptr.target_pos[i]
    inc(i)
  algorithm.sort(d_coor, system.cmp[seq_coor_t], order=algorithm.SortOrder.Ascending)
  s = 0
  e = 0
  max_s = - 1
  max_e = - 1
  max_span = - 1
  delta = seq_coor_t((0.05 * float64(max_q + max_t)))
  d_len = km_ptr.count
  d_s = - 1
  d_e = - 1
  while true:
    d_s = d_coor[s]
    d_e = d_coor[e]
    while d_e < d_s + delta and e < d_len - 1:
      inc(e, 1)
      d_e = d_coor[e]
    if max_span == - 1 or e - s > max_span:
      max_span = e - s
      max_s = s
      max_e = e
    inc(s, 1)
    if s == d_len or e == d_len:
      break
  if max_s == - 1 or max_e == - 1 or max_e - max_s < 32:
    arange.s1 = 0
    arange.e1 = 0
    arange.s2 = 0
    arange.e2 = 0
    arange.score = 0
    return arange
  newSeq(last_hit, km_ptr.count)
  newSeq(hit_score, km_ptr.count)
  newSeq(hit_count, km_ptr.count)
  i = 0
  while i < km_ptr.count:
    last_hit[i] = - 1
    hit_score[i] = 0
    hit_count[i] = 0
    inc(i)
  max_hit_idx = - 1
  max_hit_score = 0
  i = 0
  while i < km_ptr.count:
    cx = km_ptr.query_pos[i]
    cy = km_ptr.target_pos[i]
    d = cx - cy
    if d < d_coor[max_s] or d > d_coor[max_e]:
      inc(i)
      continue
    j = i - 1
    candidate_idx = - 1
    max_d = 65535
    while true:
      if j < 0: break
      px = km_ptr.query_pos[j]
      py = km_ptr.target_pos[j]
      d = px - py
      if d < d_coor[max_s] or d > d_coor[max_e]:
        dec(j)
        continue
      if cx - px > 320:
        break
      if cy > py and cx - px + cy - py < max_d and cy - py <= 320:
        max_d = cx - px + cy - py
        candidate_idx = j
      dec(j)
    if candidate_idx != - 1:
      last_hit[i] = candidate_idx
      hit_score[i] = hit_score[candidate_idx] + (64 - max_d)
      hit_count[i] = hit_count[candidate_idx] + 1
      if hit_score[i] < 0:
        hit_score[i] = 0
        hit_count[i] = 0
    else:
      hit_score[i] = 0
      hit_count[i] = 0
    if hit_score[i] > max_hit_score:
      max_hit_score = hit_score[i]
      max_hit_count = hit_count[i]
      max_hit_idx = i
    inc(i)
  if max_hit_idx == - 1:
    arange.s1 = 0
    arange.e1 = 0
    arange.s2 = 0
    arange.e2 = 0
    arange.score = 0
    return arange
  arange.score = max_hit_count + 1
  arange.e1 = km_ptr.query_pos[max_hit_idx]
  arange.e2 = km_ptr.target_pos[max_hit_idx]
  i = max_hit_idx
  while last_hit[i] != - 1:
    i = last_hit[i]
  arange.s1 = km_ptr.query_pos[i]
  arange.s2 = km_ptr.target_pos[i]
  return arange

#proc free_aln_range*(arange: ptr aln_range) =
#  discard #free(arange)
