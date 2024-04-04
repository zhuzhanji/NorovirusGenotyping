import sys
import utils as utl
import config as cg

if __name__ == "__main__":
  # pre-align GI
  utl.Clustalw(sequences = cg.paths['prefix']['GI'][0] + '.orf', output = cg.paths['prefix']['GI'][0] + '.aln')
  utl.Clustalw(sequences = cg.paths['prefix']['GI'][1] + '.orf', output = cg.paths['prefix']['GI'][1] + '.aln')
  
  # pre-align GII
  path2 = cg.paths['prefix']['GII'][0]
  utl.Clustalw(sequences = path2 + '.orf', output = path2 + '.aln')
  path2 = cg.paths['prefix']['GII'][1]
  utl.Clustalw(sequences = path2 + '.orf', output = path2 + '.aln')
  
  # pre-align all orf1, pre-align all orf2
  path2 = "./../new_data/reference_sequences_RdRp/all"
  utl.Clustalw(sequences = path2 + '.orf', output = path2 + '.aln')
  path2 = "./../new_data/reference_sequences_VP1/all"
  utl.Clustalw(sequences = path2 + '.orf', output = path2 + '.aln')
