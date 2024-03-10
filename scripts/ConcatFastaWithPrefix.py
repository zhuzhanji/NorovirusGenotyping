import sys
from Bio import SeqIO
from Bio.Seq import Seq
import json

def ConcatFastaWithPrefix(path_type, outpath: str):
    file = open(outpath, 'w')
    seq = set()
    for path, type in path_type.items():
        for seq_record in SeqIO.parse(path, "fasta"):
            name, sequence = seq_record.id, str(seq_record.seq)
            if name in seq:
              continue
            seq.add(name)
            file.write('>'+type + '##' + name+'\n')
            file.write(sequence.__str__()+'\n')
    file.close()

if __name__ == "__main__":
  contents = ''
  with open(sys.argv[1]) as f:
    contents = f.read()

  path_type = json.loads(contents)
  outpath = sys.argv[2]
  ConcatFastaWithPrefix(path_type, outpath)