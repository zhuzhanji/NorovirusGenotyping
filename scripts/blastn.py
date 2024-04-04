import sys
import utils as utl


if __name__ == "__main__":
  qrpath = sys.argv[1]
  outpath = sys.argv[2]
  print(qrpath, outpath)
  utl.QueryNcl(qrpath, outpath)

