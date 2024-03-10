import sys
from Bio import SeqIO
from Bio.Seq import Seq

def ParseFasta(inpath: str, outdir: str, option: str):
    if option == 'div':
        for seq_record in SeqIO.parse(inpath, "fasta"):
            name, sequence = seq_record.id, str(seq_record.seq)
            file = open(outdir + '/' + name + '.fasta', 'w')
            file.write(name+'\n')
            file.write(sequence+'\n')
            file.close()
    elif option == 'com':
        file = open(inpath + '.com', 'w')
        for seq_record in SeqIO.parse(inpath, "fasta"):
            name, sequence = seq_record.id, str(seq_record.seq)
            my_seq = Seq(sequence)
            my_seq_com = my_seq.complement()
            file.write('>'+name+'\n')
            file.write(my_seq_com.__str__()+'\n')
        file.close()
    elif option == 'rcom':
        file = open(inpath + '.rcom', 'w')
        for seq_record in SeqIO.parse(inpath, "fasta"):
            name, sequence = seq_record.id, str(seq_record.seq)
            my_seq = Seq(sequence)
            my_seq_com = my_seq.reverse_complement()
            file.write('>'+name+'\n')
            file.write(my_seq_com.__str__()+'\n')
        file.close()
    

if __name__ == "__main__":
  inpath = sys.argv[1]
  outdir = sys.argv[2]
  option = sys.argv[3]
  print('inpath:', inpath)
  print('outdir:', outdir)
  print('option:', option)
  ParseFasta(inpath, outdir, option)

