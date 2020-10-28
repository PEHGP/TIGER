#!/usr/bin/env python2.7
#coding=utf-8
import numpy as np
import os
import CRISPResso2
from CRISPResso2 import CRISPRessoShared
from CRISPResso2 import CRISPResso2Align
from CRISPResso2 import CRISPRessoCOREResources
def get_refobj(ref_seq,guide_seq,args):
  refs={}
  this_seq=ref_seq
  this_seq_length=len(ref_seq)
  this_name="kuan_ref"
  this_gap_incentive = np.zeros(this_seq_length+1,dtype=np.int)
  this_quant_window_coordinates = None
  guides=[guide_seq]
  (this_sgRNA_sequences, this_sgRNA_intervals, this_sgRNA_cut_points, this_sgRNA_plot_idxs, this_include_idxs,this_exclude_idxs) = CRISPRessoShared.get_amplicon_info_for_guides(this_seq,guides,args.quantification_window_center,args.quantification_window_size,this_quant_window_coordinates,args.exclude_bp_from_left,args.exclude_bp_from_right,args.plot_window_size)
  print this_sgRNA_cut_points
  for cut_point in this_sgRNA_cut_points:
    this_gap_incentive[cut_point+1] = args.needleman_wunsch_gap_incentive
  seq_rc = CRISPRessoShared.reverse_complement(this_seq)
  seeds = []
  rc_seeds = []
  seedStarts = list(range(args.exclude_bp_from_left,this_seq_length-args.exclude_bp_from_right-args.aln_seed_len,args.aln_seed_count)) #define all possible seed starts
  for seedStart in seedStarts:
      attemptsToFindSeed = 0
      thisSeedStart = seedStart
      potentialSeed = this_seq[thisSeedStart:thisSeedStart+args.aln_seed_len]
      #seed shouldn't be in reverse complement of sequence, and it should also not be the same as another seed
      while potentialSeed in seq_rc or potentialSeed in seeds:
          attemptsToFindSeed += 1
          if attemptsToFindSeed > 100:
              #raise CRISPRessoShared.AlignmentException("Can't find alignment seed that is unique to the forward sequence")
              break
          # if this seed would extend past the end of the sequence, reset to position 0 (possibly in the excluded region, but hey, we're desperate)
          if (seedStart > this_seq_length - args.aln_seed_len):
              thisSeedStart = 0
          thisSeedStart += 1
          potentialSeed = this_seq[thisSeedStart:thisSeedStart+args.aln_seed_len]
      seed_rc = CRISPRessoShared.reverse_complement(potentialSeed)
      if seed_rc in this_seq:
          #raise CRISPRessoShared.AlignmentException("Reverse compliment of seed %s is in amplicon %s even though seed is not in reverse compliment of amplicon"%(seed,amplicon))
          continue
      if potentialSeed not in seq_rc:
          seeds.append(potentialSeed)
          rc_seeds.append(seed_rc)
  this_min_aln_score = args.default_min_aln_score
  refObj = {'name':this_name,
         'sequence':this_seq,
         'sequence_length':this_seq_length,
         'min_aln_score':this_min_aln_score,
         'gap_incentive':this_gap_incentive,
         'sgRNA_cut_points':this_sgRNA_cut_points,
         'sgRNA_intervals':this_sgRNA_intervals,
         'sgRNA_sequences':this_sgRNA_sequences,
         'sgRNA_plot_idxs':this_sgRNA_plot_idxs,
         'include_idxs':this_include_idxs,
         'exclude_idxs':this_exclude_idxs,
         'idx_cloned_from':None,
         'fw_seeds':seeds,
         'rc_seeds':rc_seeds,
         }
  refs[this_name] = refObj
  return refs
def process_fastq(fastq_seq,refs,args):
    _ROOT = "/data/serverSoft/serverSoft/Anaconda/lib/python2.7/site-packages/CRISPResso2-2.0.30-py2.7-linux-x86_64.egg/CRISPResso2/"
    aln_matrix_loc = os.path.join(_ROOT,args.needleman_wunsch_aln_matrix_loc)
    CRISPRessoShared.check_file(aln_matrix_loc)
    aln_matrix = CRISPResso2Align.read_matrix(aln_matrix_loc)
    not_aln = {} #cache for reads that don't align
    count_seed_fw = 0
    count_seed_rv = 0
    count_seed_both = 0
    payload = []
    aln_scores = []
    best_match_score = -1
    best_match_s1s = []
    best_match_s2s = []
    best_match_names = []
    for idx,ref_name in enumerate(refs.keys()):
      #get alignment and score from cython
      #score = 100 * #matchedBases / length(including gaps)
      seed_i = 0
      found_forward_count = 0
      found_reverse_count = 0
      while seed_i < args.aln_seed_count and seed_i < len(refs[ref_name]['fw_seeds']):
        if refs[ref_name]['fw_seeds'][seed_i] in fastq_seq: #is forward
            found_forward_count += 1
        if refs[ref_name]['rc_seeds'][seed_i] in fastq_seq: #is rc
            found_reverse_count += 1
        seed_i += 1
      if found_forward_count > args.aln_seed_min and found_reverse_count == 0:
        fws1,fws2,fwscore=CRISPResso2Align.global_align(fastq_seq, refs[ref_name]['sequence'],matrix=aln_matrix,gap_incentive=refs[ref_name]['gap_incentive'],gap_open=args.needleman_wunsch_gap_open,gap_extend=args.needleman_wunsch_gap_extend,)
        s1 = fws1
        s2 = fws2
        score = fwscore
        count_seed_fw += 1
      elif found_forward_count == 0 and found_reverse_count > args.aln_seed_min:
        rvs1,rvs2,rvscore=CRISPResso2Align.global_align(CRISPRessoShared.reverse_complement(fastq_seq), refs[ref_name]['sequence'],matrix=aln_matrix,gap_incentive=refs[ref_name]['gap_incentive'],gap_open=args.needleman_wunsch_gap_open,gap_extend=args.needleman_wunsch_gap_extend,)
        s1 = rvs1
        s2 = rvs2
        score = rvscore
        count_seed_rv += 1
      else:
        count_seed_both += 1
        fws1,fws2,fwscore=CRISPResso2Align.global_align(fastq_seq, refs[ref_name]['sequence'],matrix=aln_matrix,gap_incentive=refs[ref_name]['gap_incentive'],gap_open=args.needleman_wunsch_gap_open,gap_extend=args.needleman_wunsch_gap_extend,)
        rvs1,rvs2,rvscore=CRISPResso2Align.global_align(CRISPRessoShared.reverse_complement(fastq_seq), refs[ref_name]['sequence'],matrix=aln_matrix,gap_incentive=refs[ref_name]['gap_incentive'],gap_open=args.needleman_wunsch_gap_open,gap_extend=args.needleman_wunsch_gap_extend,)
        s1 = fws1
        s2 = fws2
        score = fwscore
        if (rvscore > fwscore):
          s1 = rvs1
          s2 = rvs2
          score = rvscore
      aln_scores.append(score)

      #reads are matched to the reference to which they best align. The 'min_aln_score' is calculated using only the changes in 'include_idxs'
      if score > best_match_score and score > refs[ref_name]['min_aln_score']:
        best_match_score = score
        best_match_s1s = [s1]
        best_match_s2s = [s2]
        best_match_names = [ref_name]
      elif score == best_match_score:
        best_match_s1s.append(s1)
        best_match_s2s.append(s2)
        best_match_names.append(ref_name)
    if best_match_score > 0:
      #N_COMPUTED_ALN+=1
      class_names = []
      print best_match_names,best_match_s1s,best_match_s2s
      for idx in range(len(best_match_names)):
        best_match_name = best_match_names[idx]
        print refs[best_match_name]['include_idxs']
        payload=CRISPRessoCOREResources.find_indels_substitutions(best_match_s1s[idx],best_match_s2s[idx],refs[best_match_name]['include_idxs'])
        payload['ref_name'] = best_match_name
        payload['aln_scores'] = aln_scores
        is_modified = False
        if payload['deletion_n'] > 0:
          is_modified = True
        elif payload['insertion_n'] > 0:
          is_modified = True
        elif payload['substitution_n'] > 0:
          is_modified = True

        if is_modified:
          class_names.append(best_match_name+"_MODIFIED")
          payload['classification'] = 'MODIFIED'
        else:
          class_names.append(best_match_name+"_UNMODIFIED")
          payload['classification'] = 'UNMODIFIED'

        payload['aln_seq'] = best_match_s1s[idx]
        payload['aln_ref'] = best_match_s2s[idx]
        print payload
    else:
      N_COMPUTED_NOTALN+=1
      not_aln[fastq_seq] = 1
class get_args(object):
  def __init__(self,):
    super(get_args, self).__init__()
    self.needleman_wunsch_gap_incentive=1
    self.quantification_window_center=-3
    self.quantification_window_size=12
    self.exclude_bp_from_left=10
    self.exclude_bp_from_right=10
    self.plot_window_size=20
    self.needleman_wunsch_aln_matrix_loc="EDNAFULL"
    self.needleman_wunsch_gap_open=-20
    self.needleman_wunsch_gap_extend=-2
    self.aln_seed_len=10
    self.aln_seed_count=5
    self.aln_seed_min=2
    self.default_min_aln_score=60
    self.orientation=1
def main():
  ref_seq="TACACGACGCTCTTCCGATCTGTCTGTCTGTCCTCTGCCACAGGGTGCCGAGACCAGGAAGCGCTCGCCTACACTCAGCAGCCAGTTCAAGCGGTCTCTGGAGCTGCTGATGCGCACACTGGGCGCCTGCCAGCCCTTCTTTGTGCGTTGAGATCGGAAGAGCACACGTCT"
  #guide_seq="CCAGAGACCGCTTGAAC" #"GTTCAAGCGGTCTCTGG"
  guide_seq="GTTCAAGCGGTCTCTGG"
  fastq_seq="TCCCTACACGACGCTCTTCCGATCTGTCTGTCCTCTGCCACAGGGTGCCGAGACCAGGAAGCGCTCGCCTACACTCAGCAGCCAGTTCAAGCGGTCTCTGGAGCTGCTGATGCGCACACTGGGCGCCTGCCAGCCCTTCTTTGTGCGTTGAGATCGGAAGAGCACACGTCTGAAC"
  args=get_args()
  refs=get_refobj(ref_seq,guide_seq,args)
  process_fastq(fastq_seq,refs,args)
if __name__ == '__main__':
  main()
