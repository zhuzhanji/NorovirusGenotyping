import sys
import utils as utl


def FilterResults(filepath: str, outpath: str):
  results = pd.read_csv(filepath, sep="\t", header=None)
  # Define column headers
  headers = ['query', 'subject',
           'pc_identity', 'aln_length', 'mismatches', 'gaps_opened',
           'query_start', 'query_end', 'subject_start', 'subject_end',
           'e_value', 'bitscore']
  print('result:', results.shape)
  results.to_csv('out.csv', index=False)  


if __name__ == "__main__":
  qrpath = sys.argv[1]
  outpath = sys.argv[2]
  print(qrpath, outpath)
  utl.QueryNcl(qrpath, outpath)
  utl.PostQuery(outpath)

