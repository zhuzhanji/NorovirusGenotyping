import sys
import utils as utl
import config as cg

group = 'GII'
ori_seq = utl.ReadFastaFileAsDict(cg.paths['fasta'][group][0])
ori_seq.update(utl.ReadFastaFileAsDict(cg.paths['fasta'][group][1]))

viru_seq_rdrp = utl.ReadVirulignSeq(cg.paths['virulign'][group][0], skip = 1)
viru_seq_vp1 = utl.ReadVirulignSeq(cg.paths['virulign'][group][1], skip = 1)
utl.ExportAnnotation(ori_seq, viru_seq_rdrp, viru_seq_vp1, cg.paths['annotation'][group])   
