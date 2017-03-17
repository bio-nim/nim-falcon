# vim: sw=2 ts=2 sts=2 tw=80 et:
## #
## #  =====================================================================================
## # 
## #        Filename:  fastcon.c
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
import kmer_lookup_c
import DW_banded

const
  UINT8_MAX = 255
  UINT16_MAX = uint16.high.int32
  UINT_MAX = uint32.high

  MSA_BASE_GROUP_LEN = 5

type
  align_tag_t* = object
    t_pos*: seq_coor_t
    delta*: uint8
    q_base*: char
    p_t_pos*: seq_coor_t       ## # the tag position of the previous base
    p_delta*: uint8            ## # the tag delta of the previous base
    p_q_base*: char            ## # the previous base
    q_id*: uint32

  align_tags_t* = object
    length*: seq_coor_t
    align_tags*: ptr align_tag_t

  align_tag_col_t* = object
    size*: uint16
    n_link*: uint16
    #p_t_pos*: ptr seq_coor_t    ## # the tag position of the previous base
    #p_delta*: ptr uint8         ## # the tag delta of the previous base
    #p_q_base*: ptr char         ## # the previous base
    #link_count*: ptr uint16
    p_t_pos*: seq[seq_coor_t]    ## # the tag position of the previous base
    p_delta*: seq[uint8]         ## # the tag delta of the previous base
    p_q_base*: seq[char]         ## # the previous base
    link_count*: seq[uint16]
    count*: uint16
    best_p_t_pos*: seq_coor_t
    best_p_delta*: uint8
    best_p_q_base*: uint8      ## # encoded base
    score*: cdouble

  msa_base_group_t* = object
    base*: array[0 .. <MSA_BASE_GROUP_LEN, align_tag_col_t] # array of length 5

  msa_delta_group_t* = object
    # TODO(CD): size and max_delta may be redundant now
    gsize*: uint8
    #max_delta*: uint8
    delta*: seq[msa_base_group_t]

  msa_pos_t* = seq[msa_delta_group_t]

common.usePtr[align_tag_t]
common.usePtr[ptr align_tags_t]
common.usePtr[align_tag_col_t]
common.usePtr[msa_pos_t]
common.usePtr[msa_base_group_t]
common.usePtr[seq_coor_t]
common.usePtr[uint16]
common.usePtr[uint8]
common.usePtr[int8]
common.usePtr[char]


proc get_align_tags*(aln_q_seq: cstring; aln_t_seq: cstring; aln_seq_len: seq_coor_t;
                    srange: ref aln_range; q_id: uint32; t_offset: seq_coor_t): ref align_tags_t =
  var p_q_base: char
  var tags: ref align_tags_t
  var
    i: seq_coor_t
    j: seq_coor_t
    jj: seq_coor_t
    k: seq_coor_t
    p_j: seq_coor_t
    p_jj: seq_coor_t
  new(tags)
  tags.length = aln_seq_len
  tags.align_tags = calloc[align_tag_t](aln_seq_len + 1)
  i = srange.s1 - 1
  j = srange.s2 - 1
  jj = 0
  p_j = - 1
  p_jj = 0
  p_q_base = '.'
  k = 0
  while k < aln_seq_len:
    if aln_q_seq[k] != '-':
      inc(i)
      inc(jj)
    if aln_t_seq[k] != '-':
      inc(j)
      jj = 0
    if j + t_offset >= 0 and jj < UINT8_MAX and p_jj < UINT8_MAX:
      (tags.align_tags[k]).t_pos = j + t_offset
      (tags.align_tags[k]).delta = uint8(jj)
      (tags.align_tags[k]).p_t_pos = p_j + t_offset
      (tags.align_tags[k]).p_delta = uint8(p_jj)
      (tags.align_tags[k]).p_q_base = p_q_base
      (tags.align_tags[k]).q_base = aln_q_seq[k]
      (tags.align_tags[k]).q_id = q_id
      p_j = j
      p_jj = jj
      p_q_base = aln_q_seq[k]
    inc(k)
  ## # sentinal at the end
  ## #k = aln_seq_len;
  tags.length = k
  (tags.align_tags[k]).t_pos = -1 #UINT_MAX
  (tags.align_tags[k]).delta = UINT8_MAX
  (tags.align_tags[k]).q_base = '.'
  (tags.align_tags[k]).q_id = UINT_MAX
  return tags

proc free_align_tags*(tags: ref align_tags_t) =
  dealloc(tags.align_tags)
  #dealloc(tags)

proc allocate_aln_col*(col: ptr align_tag_col_t) =
  newSeq(col.p_t_pos, col.size)
  newSeq(col.p_delta, col.size)
  newSeq(col.p_q_base, col.size)
  newSeq(col.link_count, col.size)
  #col.p_t_pos = calloc0[seq_coor_t](col.size)
  #col.p_delta = calloc0[uint8](col.size)
  #col.p_q_base = calloc0[char](col.size)
  #col.link_count = calloc0[uint16](col.size)

proc realloc_aln_col*(col: ptr align_tag_col_t) =
  #echo "realloc_aln_col @", cast[ByteAddress](col), " to ", col.size
  col.p_t_pos.add(newSeq[seq_coor_t](col.size))
  col.p_delta.add(newSeq[uint8](col.size))
  col.p_q_base.add(newSeq[char](col.size))
  col.link_count.add(newSeq[uint16](col.size))
  #col.p_t_pos = realloc0[seq_coor_t](col.p_t_pos, col.size)
  #col.p_delta = realloc0[uint8](col.p_delta, col.size)
  #col.p_q_base = realloc0[char](col.p_q_base, col.size)
  #col.link_count = realloc0[uint16](col.link_count, col.size)

proc free_aln_col*(col: ptr align_tag_col_t) =
  discard
  #dealloc(col.p_t_pos)
  #dealloc(col.p_delta)
  #dealloc(col.p_q_base)
  #dealloc(col.link_count)

proc allocate_delta_group*(g: ptr msa_delta_group_t, cap: int) =
  var
    i: uint8
    j: cint
  #g.max_delta = 0 # -1?
  #g.delta = calloc0[msa_base_group_t](g.size)
  g.delta = newSeqOfCap[msa_base_group_t](cap)
  i = 0
  g.gsize = 0
  while i < g.gsize:
    #g.delta[i].base = calloc[align_tag_col_t](5)
    j = 0
    while j < 5:
      g.delta[i.int].base[j].size = 8
      allocate_aln_col(addr((g.delta[i.int].base[j])))
      inc(j)
    inc(i)

proc realloc_delta_group*(g: ptr msa_delta_group_t; new_size: int) =
  #log("realloc_del:", $new_size, " was:", $g.delta.len)
  assert(new_size <= 255)
  # new_size was uint8
  var
    i: int #uint16
    j: int #cint
    bs: int #uint8
    es: int #uint16
  bs = g.gsize.int
  es = new_size
  g.delta.setLen(new_size) # calloc(..., new_size)
  i = bs
  while i < es:
    #g.delta[i].base = calloc[align_tag_col_t](5)
    j = 0
    while j < 5:
      g.delta[i].base[j].size = 8
      allocate_aln_col(addr((g.delta[i].base[j])))
      inc(j)
    inc(i)
  g.gsize= new_size.uint8

proc free_delta_group*(g: ptr msa_delta_group_t) =
  discard """
  var
    i: uint8
    j: cint
  i = 0
  while i < g.size:
    j = 0
    while j < 5:
      free_aln_col(addr(g.delta[i].base[j]))
      inc(j)
    #dealloc(g.delta[i].base)
    inc(i)
  dealloc(g.delta)
  """

proc update_col*(col: ptr align_tag_col_t; p_t_pos: seq_coor_t; p_delta: uint8;
                p_q_base: char) =
  var updated: cint = 0
  var kk: int16
  inc(col.count, 1)
  kk = 0
  while kk < col.n_link.int16:
    if p_t_pos == col.p_t_pos[kk.int] and p_delta == col.p_delta[kk.int] and
        p_q_base == col.p_q_base[kk]:
      inc(col.link_count[kk])
      updated = 1
      break
    inc(kk)
  if updated == 0:
    #log("col.n_link=", $col.n_link)
    if col.n_link + 1 > col.size:
      #log("1col.size=", $col.size)
      if col.size < uint16((UINT16_MAX shr 1) - 1):
        col.size = col.size * 2
      else:
        inc(col.size, 256)
      assert(col.size < uint16(UINT16_MAX - 1))
      realloc_aln_col(col)
    #log("repr(col)=", repr(col))
    kk = col.n_link.int16
    #log("kk=", $kk)
    col.p_t_pos[kk.int] = p_t_pos
    col.p_delta[kk] = p_delta
    col.p_q_base[kk] = p_q_base
    col.link_count[kk] = 1
    inc(col.n_link)

proc get_msa_working_sapce*(max_t_len: cuint): ref msa_pos_t =
  var msa_array: ref msa_pos_t
  new(msa_array)
  newSeq(msa_array[], max_t_len)
  var i: int
  #msa_array = calloc[msa_pos_t](max_t_len)
  i = 0
  while i < max_t_len.int:
    #msa_array[i] = calloc[msa_delta_group_t](1)
    msa_array[i].gsize = 8
    allocate_delta_group(addr msa_array[i], 8)
    inc(i)
  return msa_array

proc clean_msa_working_space*(msa_array: ref msa_pos_t; max_t_len: int) =
  # max_t_len was uint
  var
    i: int #uint
    j: int #uint
    k: int #uint
  i = 0
  while i < max_t_len:
    realloc_delta_group(addr msa_array[i], 0)
    discard """
    j = 0
    while j < msa_array[i].delta.len
      k = 0
      while k < 5:
        let col: ptr align_tag_col_t = addr msa_array[i].delta[j].base[k]
        ## #
        ## #                for (c =0; c < col->size; c++) {
        ## #                    col->p_t_pos[c] = 0;
        ## #                    col->p_delta[c] = 0;
        ## #                    col->p_q_base[c] = 0;
        ## #                    col->link_count[c] =0;
        ## #                }
        ## #                
        col.n_link = 0
        col.count = 0
        col.best_p_t_pos = 0
        col.best_p_delta = 0
        col.best_p_q_base = 0
        col.score = 0
        inc(k)
      inc(j)
    msa_array[i].max_delta = 0
      """
    inc(i)

#const
#  STATIC_ALLOCATE* = true

when defined(STATIC_ALLOCATE):
  log("STATIC_ALLOCATE")
  var msa_array {.threadvar.}: ref msa_pos_t


proc get_cns_from_align_tags*(tag_seqs: seq[ref align_tags_t]; n_tag_seqs: seq_coor_t;
                             t_len: seq_coor_t; min_cov: int): ref consensus_data =
  var
    i: seq_coor_t
    j: seq_coor_t
  var t_pos: seq_coor_t = 0
  var coverage: seq[int]
  var local_nbase: seq[uint8]
  var consensus: ref consensus_data
  ## #char * consensus;
  var c_tag: ptr align_tag_t
  newSeq(coverage, t_len)
  newSeq(local_nbase, t_len)
  when not defined(STATIC_ALLOCATE):
    log("no static_ALLOCATE")
    var msa_array: ref msa_pos_t
    new(msa_array)
    #msa_array = calloc[msa_pos_t](t_len)
    newSeq(msa_array[], t_len)
    i = 0
    while i < t_len:
      #msa_array[i] = calloc[msa_delta_group_t](1)
      msa_array[i].size = 8
      allocate_delta_group(addr msa_array[i], 8)
      inc(i)
  else:
    assert(t_len < 100000)
    if msa_array == nil:
      msa_array = get_msa_working_sapce(100000)
  ## # loop through every alignment
  i = 0
  while i < n_tag_seqs:
    ## # for each alignment position, insert the alignment tag to msa_array
    j = 0
    while j < tag_seqs[i].length:
      c_tag = tag_seqs[i].align_tags + j
      var delta: uint8
      delta = c_tag.delta
      if delta == 0:
        t_pos = c_tag.t_pos
        inc(coverage[t_pos])
      if msa_array[t_pos].delta.len <= delta.int:
        #msa_array[t_pos].max_delta = delta
        #if msa_array[t_pos].max_delta > msa_array[t_pos].delta.len.uint8:
        realloc_delta_group(addr msa_array[t_pos], delta.int + 1)
      var base: uint8
      case c_tag.q_base
      of 'A':
        base = 0
      of 'C':
        base = 1
      of 'G':
        base = 2
      of 'T':
        base = 3
      of '-':
        base = 4
      else:
        base = 255
      ## # Note: On bad input, base may be -1.
      #log("tpos:", $tpos, " len(arr):", $len(msa_array[]))
      #log("delta:", $delta, " len(arr.delta]):", $msa_array[t_pos].delta.len)
      #log("repr(msa_array[t_pos].delta[delta].base)@", $base, "=", $repr(msa_array[t_pos].delta[delta.int].base))
      update_col(addr((msa_array[t_pos].delta[delta.int].base[base])), c_tag.p_t_pos,
                 c_tag.p_delta, c_tag.p_q_base)
      inc(local_nbase[t_pos])
      inc(j)
    inc(i)
  ## # propogate score throught the alignment links, setup backtracking information
  var g_best_aln_col: ptr align_tag_col_t = nil
  var g_best_ck: int16 = 0
  var g_best_t_pos: seq_coor_t = 0
  var kk: cint
  ## # char base;
  var g_best_score: cdouble
  ## # char best_mark;
  g_best_score = - 1
  i = 0
  while i < t_len:
    ## #loop through every template base
    ## #printf("max delta: %d %d\n", i, msa_array[i]->max_delta);
    j = 0
    while j < seq_coor_t(msa_array[i].delta.len):
      ## # loop through every delta position
      kk = 0
      while kk < MSA_BASE_GROUP_LEN:
        ## # loop through diff bases of the same delta posiiton
        ## #
        ## #                    switch (kk) {
        ## #                        case 0: base = 'A'; break;
        ## #                        case 1: base = 'C'; break;
        ## #                        case 2: base = 'G'; break;
        ## #                        case 3: base = 'T'; break;
        ## #                        case 4: base = '-'; break;
        ## #                    }
        ## #                    
        let aln_col: ptr align_tag_col_t = addr msa_array[i].delta[j].base[kk]
        if aln_col.count >= 0'u16:
          var score: cdouble
          var best_score: cdouble = -1
          var best_i: int = -1
          var best_j: int = 255 #uint8
          var best_b: int = 255 #uint8
          var best_ck: int16 = 255 #int16
          var ck: int16 = 0
          while ck < aln_col.n_link.int16:
            ## # loop through differnt link to previous column
            var pi: int #cint
            var pj: int #uint8
            var pkk: int #uint8
            pi = aln_col.p_t_pos[ck.int]
            pj = aln_col.p_delta[ck].int
            case aln_col.p_q_base[ck]
            of 'A':
              pkk = 0
            of 'C':
              pkk = 1
            of 'G':
              pkk = 2
            of 'T':
              pkk = 3
            of '-':
              pkk = 4
            else:
              pkk = 4
            if aln_col.p_t_pos[ck.int] == - 1:
              score = cdouble(aln_col.link_count[ck]) -
                  cdouble(coverage[i]) * 0.5
            else:
              score = msa_array[pi].delta[pj].base[pkk].score +
                  cdouble(aln_col.link_count[ck]) -
                  cdouble(coverage[i]) * 0.5
            ## # best_mark = ' ';
            if score > best_score:
              best_score = score
              best_i = pi
              aln_col.best_p_t_pos = best_i.seq_coor_t
              best_j = pj
              aln_col.best_p_delta = best_j.uint8
              best_b = pkk
              aln_col.best_p_q_base = best_b.uint8
              best_ck = ck
              ## # best_mark = '*';
            inc(ck)
          aln_col.score = best_score
          if best_score > g_best_score:
            g_best_score = best_score
            g_best_aln_col = aln_col
            g_best_ck = best_ck
            g_best_t_pos = i
            ## #printf("GB %d %d %d %d\n", i, j, ck, g_best_aln_col);
        inc(kk)
      inc(j)
    inc(i)
  assert(g_best_score != - 1)
  ## # reconstruct the sequences
  var index: seq_coor_t
  var bb: char = '$'
  var ck: range[0 .. 4]
  var cns_str: ptr char #seq[char]
  var cns_int: ptr int8 #seq[int8]
  var eqv: ptr seq[cint]
  var score0: cdouble
  new(consensus)
  #newSeq(consensus.sequence, t_len*2+1)
  consensus.sequence = newString(t_len*2+1)
  newSeq(consensus.eqv, tlen*2+1)
  # probably do not need +1

  cns_str = addr consensus.sequence[0] # alias
  cns_int = cast[ptr int8](addr consensus.sequence[0]) # alias
  eqv = addr consensus.eqv # alias
  index = 0
  ck = g_best_ck
  i = g_best_t_pos
  while true:
    if coverage[i] > min_cov:
      case ck
      of 0:
        bb = 'A'
      of 1:
        bb = 'C'
      of 2:
        bb = 'G'
      of 3:
        bb = 'T'
      of 4:
        bb = '-'
    else:
      case ck
      of 0:
        bb = 'a'
      of 1:
        bb = 'c'
      of 2:
        bb = 'g'
      of 3:
        bb = 't'
      of 4:
        bb = '-'
    ## # Note: On bad input, bb will keep previous value, possibly '$'.
    score0 = g_best_aln_col.score
    i = g_best_aln_col.best_p_t_pos
    if i == - 1 or index >= t_len * 2: break
    j = seq_coor_t(g_best_aln_col.best_p_delta)
    ck = g_best_aln_col.best_p_q_base
    g_best_aln_col = addr msa_array[i].delta[j].base[ck]
    if bb != '-':
      cns_str[index] = bb
      eqv[index] = cint(score0) - cint(g_best_aln_col.score)
      ## #printf("C %d %d %c %lf %d %d\n", i, index, bb, g_best_aln_col->score, coverage[i], eqv[index] );
      inc(index)
  ## # reverse the sequence
  i = 0
  while i < index div 2:
    cns_int[i] = cns_int[i] xor cns_int[index - i - 1]
    cns_int[index - i - 1] = cns_int[i] xor cns_int[index - i - 1]
    cns_int[i] = cns_int[i] xor cns_int[index - i - 1]
    eqv[i] = eqv[i] xor eqv[index - i - 1]
    eqv[index - i - 1] = eqv[i] xor eqv[index - i - 1]
    eqv[i] = eqv[i] xor eqv[index - i - 1]
    inc(i)
  cns_int[index] = 0
  ## #printf("%s\n", cns_str);
  when not defined(STATIC_ALLOCATE):
    i = 0
    while i < t_len:
      free_delta_group(msa_array[i])
      dealloc(msa_array[i])
      inc(i)
    dealloc(msa_array)
  else:
    clean_msa_working_space(msa_array, (t_len + 1))
  return consensus

## #const unsigned int K = 8;

proc generate_consensus*(input_seq: cStringArray; n_seq: int; min_cov: int;
                        K: int; min_idt: float32): consensus_data =
  var j: int
  var seq_count: int
  var aligned_seq_count: seq_coor_t
  var sa_ptr: seq_array
  var sda_ptr: seq_addr_array
  var kmer_match_ptr: ref kmer_match
  var arange: ref aln_range
  var aln: ref alignment
  var tags_list: seq[ref align_tags_t]
  ## #char * consensus;
  var consensus: ref consensus_data
  var max_diff: cdouble
  max_diff = 1.0 - min_idt
  seq_count = n_seq
  ## #for (j=0; j < seq_count; j++) {
  ## #    printf("seq_len: %u %u\n", j, strlen(input_seq[j]));
  ## #};
  ## #fflush(stdout);
  #var ii: int = 0
  #while ii < seq_count:
  #  echo "seq_len:", ii, " ", len(input_seq[ii])
  #  inc(ii)
  newSeq(tags_list, seq_count)
  var lk_ptr = kmer_lookup_c.allocate_kmer_lookup(seq_coor_t(1 shl (K * 2)))
  #log("len(lk_ptr):", $(1 shl (K * 2)), "==", $len(lk_ptr))
  sa_ptr = allocate_seq(len(input_seq[0]).seq_coor_t)
  sda_ptr = allocate_seq_addr(len(input_seq[0]).seq_coor_t)
  # TODO(CD): That was 0s for .c, but not 0s here. Where should we init?
  add_sequence(0.seq_coor_t, K.cuint, input_seq[0], len(input_seq[0]).seq_coor_t, sda_ptr, sa_ptr, lk_ptr)
  ## #mask_k_mer(1 << (K * 2), lk_ptr, 16);
  aligned_seq_count = 0
  j = 1
  while j < seq_count:
    ## #printf("seq_len: %ld %u\n", j, strlen(input_seq[j]));
    kmer_match_ptr = find_kmer_pos_for_seq(input_seq[j.int], len(input_seq[j.int]).seq_coor_t, K.cuint,
        sda_ptr, lk_ptr)
    const
      INDEL_ALLOWENCE_0 = 6
    arange = find_best_aln_range(kmer_match_ptr, K.int * INDEL_ALLOWENCE_0, 5.seq_coor_t)
    #log("ARANGE:", repr(arange))
    ## # narrow band to avoid aligning through big indels
    ## #printf("1:%ld %ld %ld %ld\n", arange_->s1, arange_->e1, arange_->s2, arange_->e2);
    ## #arange = find_best_aln_range2(kmer_match_ptr, K * INDEL_ALLOWENCE_0, 5.seq_coor_t);  // narrow band to avoid aligning through big indels
    ## #printf("2:%ld %ld %ld %ld\n\n", arange->s1, arange->e1, arange->s2, arange->e2);
    const
      INDEL_ALLOWENCE_1 = 0.1
    if arange.e1 - arange.s1 < 100 or arange.e2 - arange.s2 < 100 or
        abs((arange.e1 - arange.s1) - (arange.e2 - arange.s2)) >
        (int)(0.5 * INDEL_ALLOWENCE_1 *
        float(arange.e1 - arange.s1 + arange.e2 - arange.s2)):
      inc(j)
      continue
    const
      INDEL_ALLOWENCE_2 = 150
    aln = DW_banded.align((addr input_seq[j.int][0]) + arange.s1, arange.e1 - arange.s1,
              (addr input_seq[0][0]) + arange.s2, arange.e2 - arange.s2, INDEL_ALLOWENCE_2, true)
    if (aln.aln_str_size > 500) and ((cdouble(aln.dist) / cdouble(aln.aln_str_size)) < max_diff):
      tags_list[aligned_seq_count] = get_align_tags(aln.q_aln_str, aln.t_aln_str,
          aln.aln_str_size, arange, j.uint32, 0.seq_coor_t)
      inc(aligned_seq_count)
    #free_alignment(aln)
    inc(j)
  if aligned_seq_count > 0:
    consensus = get_cns_from_align_tags(tags_list, aligned_seq_count,
                                      len(input_seq[0]).seq_coor_t, min_cov)
  else:
    ## # allocate an empty consensus sequence
    new(consensus)
    #newSeq(consensus.sequence, 1)
    consensus.sequence = newString(1)
    newSeq(consensus.eqv, 1)
  ## #free(consensus);
  j = 0
  while j < aligned_seq_count:
    free_align_tags(tags_list[j])
    inc(j)
  #free(tags_list)
  return consensus[]

proc generate_utg_consensus*(input_seq: cStringArray; in_offset: seq[seq_coor_t];
                            n_seq: int; min_cov: int; K: int; min_idt: cdouble): consensus_data =
  var offset: seq[seq_coor_t] = in_offset
  var j: int
  var seq_count: int
  var aligned_seq_count: int
  var arange: ref aln_range
  var aln: ref alignment
  var tags_list: seq[ref align_tags_t]
  ## #char * consensus;
  var consensus: ref consensus_data
  var max_diff: cdouble
  var utg_len: seq_coor_t
  var r_len: seq_coor_t
  new(consensus)
  #newSeq(consensus.sequence, 1)
  #newSeq(consensus.eqv, 1)
  max_diff = 1.0 - min_idt
  seq_count = n_seq
  ## #**
  ## #    for (j=0; j < seq_count; j++) {
  ## #        printf("seq_len: %u %u\n", j, strlen(input_seq[j]));
  ## #    };
  ## #    fflush(stdout);
  ## #    *
  newSeq(tags_list, seq_count + 1)
  utg_len = len(input_seq[0]).seq_coor_t
  aligned_seq_count = 0
  new(arange)
  arange.s1 = 0
  arange.e1 = len(input_seq[0]).seq_coor_t
  arange.s2 = 0
  arange.e2 = len(input_seq[0]).seq_coor_t
  tags_list[aligned_seq_count.int] = get_align_tags(input_seq[0], input_seq[0],
      len(input_seq[0]).seq_coor_t, arange, 0, 0)
  inc(aligned_seq_count, 1)
  j = 1
  while j < seq_count:
    arange.s1 = 0
    arange.e1 = seq_coor_t(len(input_seq[j]) - 1)
    arange.s2 = 0
    arange.e2 = seq_coor_t(len(input_seq[j]) - 1)
    r_len = len(input_seq[j]).seq_coor_t
    ## #printf("seq_len: %u %u\n", j, r_len);
    if offset[j.int] < 0:
      if (r_len + offset[j]) < 128:
        inc(j)
        continue
      if r_len + offset[j] < utg_len:
        ## #printf("1: %ld %u %u\n", offset[j], r_len, utg_len);
        aln = DW_banded.align((addr input_seq[j][0]) - offset[j], r_len + offset[j], (addr input_seq[0][0]),
                  r_len + offset[j], 500, true)
      else:
        ## #printf("2: %ld %u %u\n", offset[j], r_len, utg_len);
        aln = DW_banded.align((addr input_seq[j][0]) - offset[j], utg_len, (addr input_seq[0][0]), utg_len, 500, true)
      offset[j] = 0
    else:
      if offset[j] > utg_len - 128:
        inc(j)
        continue
      if offset[j] + r_len > utg_len:
        ## #printf("3: %ld %u %u\n", offset[j], r_len, utg_len);
        aln = DW_banded.align(addr input_seq[j][0], utg_len - offset[j], (addr input_seq[0][0]) + offset[j],
                  utg_len - offset[j], 500, true)
      else:
        ## #printf("4: %ld %u %u\n", offset[j], r_len, utg_len);
        aln = DW_banded.align(addr input_seq[j][0], r_len, (addr input_seq[0][0]) + offset[j], r_len, 500, true)
    if aln.aln_str_size > 500 and
        cdouble(aln.dist) / cdouble(aln.aln_str_size) < max_diff:
      tags_list[aligned_seq_count] = get_align_tags(aln.q_aln_str, aln.t_aln_str,
          aln.aln_str_size, arange, j.uint32, offset[j])
      inc(aligned_seq_count)
    #free_alignment(aln)
    inc(j)
  if aligned_seq_count > 0:
    consensus = get_cns_from_align_tags(tags_list, aligned_seq_count.seq_coor_t, utg_len, 0)
  else:
    ## # allocate an empty consensus sequence
    consensus.sequence = ""
    consensus.eqv = @[]
  ## #free(consensus);
  j = 0
  while j < aligned_seq_count:
    free_align_tags(tags_list[j])
    inc(j)
  #free(tags_list)
  return consensus[]

#proc free_consensus_data*(consensus: ptr consensus_data) =
#  free(consensus.sequence)
#  free(consensus.eqv)
#  free(consensus)
