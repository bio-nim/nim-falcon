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

type
  seq_coor_t* = cint
  alignment* {.importc: "alignment", header: "common.h".} = object
    aln_str_size* {.importc: "aln_str_size".}: seq_coor_t
    dist* {.importc: "dist".}: seq_coor_t
    aln_q_s* {.importc: "aln_q_s".}: seq_coor_t
    aln_q_e* {.importc: "aln_q_e".}: seq_coor_t
    aln_t_s* {.importc: "aln_t_s".}: seq_coor_t
    aln_t_e* {.importc: "aln_t_e".}: seq_coor_t
    q_aln_str* {.importc: "q_aln_str".}: cstring
    t_aln_str* {.importc: "t_aln_str".}: cstring

  d_path_data* {.importc: "d_path_data", header: "common.h".} = object
    pre_k* {.importc: "pre_k".}: seq_coor_t
    x1* {.importc: "x1".}: seq_coor_t
    y1* {.importc: "y1".}: seq_coor_t
    x2* {.importc: "x2".}: seq_coor_t
    y2* {.importc: "y2".}: seq_coor_t

  d_path_data2* {.importc: "d_path_data2", header: "common.h".} = object
    d* {.importc: "d".}: seq_coor_t
    k* {.importc: "k".}: seq_coor_t
    pre_k* {.importc: "pre_k".}: seq_coor_t
    x1* {.importc: "x1".}: seq_coor_t
    y1* {.importc: "y1".}: seq_coor_t
    x2* {.importc: "x2".}: seq_coor_t
    y2* {.importc: "y2".}: seq_coor_t

  path_point* {.importc: "path_point", header: "common.h".} = object
    x* {.importc: "x".}: seq_coor_t
    y* {.importc: "y".}: seq_coor_t

  kmer_lookup* {.importc: "kmer_lookup", header: "common.h".} = object
    start* {.importc: "start".}: seq_coor_t
    last* {.importc: "last".}: seq_coor_t
    count* {.importc: "count".}: seq_coor_t

  base* = cuchar
  seq_array* = ptr base
  seq_addr* = seq_coor_t
  seq_addr_array* = ptr seq_addr
  kmer_match* {.importc: "kmer_match", header: "common.h".} = object
    count* {.importc: "count".}: seq_coor_t
    query_pos* {.importc: "query_pos".}: ptr seq_coor_t
    target_pos* {.importc: "target_pos".}: ptr seq_coor_t

  aln_range* {.importc: "aln_range", header: "common.h".} = object
    s1* {.importc: "s1".}: seq_coor_t
    e1* {.importc: "e1".}: seq_coor_t
    s2* {.importc: "s2".}: seq_coor_t
    e2* {.importc: "e2".}: seq_coor_t
    score* {.importc: "score".}: clong

  consensus_data* {.importc: "consensus_data", header: "common.h".} = object
    sequence* {.importc: "sequence".}: cstring
    eqv* {.importc: "eqv".}: ptr cint


proc allocate_kmer_lookup*(a2: seq_coor_t): ptr kmer_lookup {.cdecl,
    importc: "allocate_kmer_lookup", header: "common.h".}
proc init_kmer_lookup*(a2: ptr kmer_lookup; a3: seq_coor_t) {.cdecl,
    importc: "init_kmer_lookup", header: "common.h".}
proc free_kmer_lookup*(a2: ptr kmer_lookup) {.cdecl, importc: "free_kmer_lookup",
    header: "common.h".}
proc allocate_seq*(a2: seq_coor_t): seq_array {.cdecl, importc: "allocate_seq",
    header: "common.h".}
proc init_seq_array*(a2: seq_array; a3: seq_coor_t) {.cdecl,
    importc: "init_seq_array", header: "common.h".}
proc free_seq_array*(a2: seq_array) {.cdecl, importc: "free_seq_array",
                                   header: "common.h".}
proc allocate_seq_addr*(size: seq_coor_t): seq_addr_array {.cdecl,
    importc: "allocate_seq_addr", header: "common.h".}
proc free_seq_addr_array*(a2: seq_addr_array) {.cdecl,
    importc: "free_seq_addr_array", header: "common.h".}
proc find_best_aln_range*(a2: ptr kmer_match; a3: seq_coor_t; a4: seq_coor_t;
                         a5: seq_coor_t): ptr aln_range {.cdecl,
    importc: "find_best_aln_range", header: "common.h".}
proc free_aln_range*(a2: ptr aln_range) {.cdecl, importc: "free_aln_range",
                                      header: "common.h".}
proc find_kmer_pos_for_seq*(a2: cstring; a3: seq_coor_t; K: cuint; a5: seq_addr_array;
                           a6: ptr kmer_lookup): ptr kmer_match {.cdecl,
    importc: "find_kmer_pos_for_seq", header: "common.h".}
proc free_kmer_match*(`ptr`: ptr kmer_match) {.cdecl, importc: "free_kmer_match",
    header: "common.h".}
proc add_sequence*(a2: seq_coor_t; a3: cuint; a4: cstring; a5: seq_coor_t;
                  a6: seq_addr_array; a7: seq_array; a8: ptr kmer_lookup) {.cdecl,
    importc: "add_sequence", header: "common.h".}
proc mask_k_mer*(a2: seq_coor_t; a3: ptr kmer_lookup; a4: seq_coor_t) {.cdecl,
    importc: "mask_k_mer", header: "common.h".}
proc align*(a2: cstring; a3: seq_coor_t; a4: cstring; a5: seq_coor_t; a6: seq_coor_t;
           a7: cint): ptr alignment {.cdecl, importc: "align", header: "common.h".}
proc free_alignment*(a2: ptr alignment) {.cdecl, importc: "free_alignment",
                                      header: "common.h".}
proc free_consensus_data*(a2: ptr consensus_data) {.cdecl,
    importc: "free_consensus_data", header: "common.h".}
proc generate_consensus*(input_seq: cstringArray; n_seq: cuint; min_cov: cuint;
                        K: cuint; min_idt: cdouble): ptr consensus_data {.cdecl,
    importc: "generate_consensus", header: "common.h".}
proc find_best_aln_range2*(km_ptr: ptr kmer_match; K: seq_coor_t;
                          bin_width: seq_coor_t; count_th: seq_coor_t): ptr aln_range {.
    cdecl, importc: "find_best_aln_range2", header: "common.h".}
