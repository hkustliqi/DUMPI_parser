#!/usr/bin/env python3

import argparse
import os
import gzip
import re
import copy
import string


#This script converts n DUMPI logs(.bin) to readable format (.txt)
def main(args):
  n = int(args.n)

  for i in range(n):
    # find all log names
    s = list(args.DUMPI)
    num = list(str(i).zfill(4))
    for j in range(4):
      s[len(s)-8+j] = num[j]
    inputName = ''.join(s)
    outputTxt = inputName.replace('bin', 'txt')
    # dumpi2ascii converts bin to txt
    os.system('LD_LIBRARY_PATH=/home/qi/sst-dumpi/lib/ /home/qi/sst-dumpi/bin/dumpi2ascii ' + inputName + ' > ' + outputTxt)
    

if __name__ == '__main__':
  ap = argparse.ArgumentParser()
  ap.add_argument('DUMPI',
                  help='any one of the DUMPI log file series (format: ****_000*.txt)')
  ap.add_argument('n',
                  help='totol number of files')
  args = ap.parse_args()
  main(args)

