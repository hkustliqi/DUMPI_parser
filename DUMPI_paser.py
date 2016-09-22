#!/usr/bin/env python3

import argparse
import os
import gzip
import re
import copy
import math
import csv

flit_size = 1 #bytes
header_size = 1 #flit
ack_size = 1 #flit

#TODO: add all MPI data types
#unit is byte
def MPI_Data_Type_to_size(datatype):
  if datatype == 2: # MPI char
    return 1
  elif datatype == 3: # MPI signed char
    return 1
  elif datatype == 4: # MPI unsigned char
    return 1
  elif datatype == 5: # MPI byte
    return 1
  elif datatype == 6: # MPI WCHAR
    return 4
  elif datatype == 7: # MPI short
    return 2
  elif datatype == 8: # MPI unsigned short
    return 2
  elif datatype == 9: # MPI int
    return 2
  elif datatype == 10: # MPI unsigned int
    return 2
  elif datatype == 11: # MPI long
    return 4
  elif datatype == 12: # MPI unsigned long
    return 4
  elif datatype == 13: # MPI float
    return 4
  elif datatype == 14: # MPI double
    return 8
  elif datatype == 15: # MPI long double
    return 10
  elif datatype == 17: # MPI unsigned long long
    return 8
  elif datatype == 18: # MPI long long
    return 8
  elif datatype == 19: # MPI packed
    return 8
  elif datatype == 28: # user defined datatype
    return 8
  else:
    print('undefined MPI datatype', datatype)

def analyze_DUMPI(fd, matrix):
  lines = fd.readlines()
  for i in range(0, len(lines)):
    line = lines[i]
    lineStr = line.decode('utf-8')

    # find MPI rank
    if 'MPI_Comm_rank entering' in lineStr:
      rankLine = lines[i+2].decode('utf-8')
      rank = re.search('rank=(\d+)', rankLine).group(1)
      src = int(rank)
    # map MPI functions to traffic
    # MPI_send and variations
    elif ('MPI_Send entering' or 'MPI_Isend entering' or 'MPI_Ssend entering' or 'MPI_Issend entering' or 'MPI_Sendrecv entering') in lineStr:
      countLine = lines[i+1].decode('utf-8')
      count = re.search('count=(\d+)', countLine).group(1)
      datatypeLine = lines[i+2].decode('utf-8')
      datatype = re.search('datatype=(\d+)', datatypeLine).group(1)
      size = MPI_Data_Type_to_size(int(datatype))
      destLine = lines[i+3].decode('utf-8')
      dest = int(re.search('dest=(\d+)', destLine).group(1))
      matrix[src][dest] += header_size + math.ceil (int(count) * size / flit_size)
      matrix[dest][src] += ack_size
    # MPI_Bcast
    elif 'MPI_Bcast entering' in lineStr:
      countLine = lines[i+1].decode('utf-8')
      count = re.search('count=(\d+)', countLine).group(1)
      datatypeLine = lines[i+2].decode('utf-8')
      datatype = re.search('datatype=(\d+)', datatypeLine).group(1)
      size = MPI_Data_Type_to_size(int(datatype))
      rootLine = lines[i+3].decode('utf-8')
      root = re.search('root=(\d+)', rootLine).group(1)
      for i in range(len(matrix[0])):
        matrix[int(root)][i] += header_size + math.ceil (int(count) * size / len(matrix) /flit_size)
        matrix[i][int(root)] += ack_size
    # MPI_Alltoall
    elif ('MPI_Alltoall entering' or 'MPI_Alltoallv entering') in lineStr:
      countLine = lines[i+1].decode('utf-8')
      count = re.search('sendcount=(\d+)', countLine).group(1)
      datatypeLine = lines[i+2].decode('utf-8')
      datatype = re.search('sendtype=(\d+)', datatypeLine).group(1)
      size = MPI_Data_Type_to_size(int(datatype))
      for i in range(len(matrix)):
        for j in range(len(matrix[0])):
          matrix[i][j] += header_size + math.ceil (int(count) * size / len(matrix) / flit_size)
          matrix[j][i] += ack_size
    # MPI_Gather
    elif ('MPI_Gather entering' or 'MPI_Gatherv entering') in lineStr:
      countLine = lines[i+2].decode('utf-8')
      count = re.search('sendcount=(\d+)', countLine).group(1)
      datatypeLine = lines[i+3].decode('utf-8')
      datatype = re.search('sendtype=(\d+)', datatypeLine).group(1)
      size = MPI_Data_Type_to_size(int(datatype))
      rootLine = lines[i+6].decode('utf-8')
      root = re.search('root=(\d+)', countLine).group(1)
      for i in range(len(matrix)):
        matrix[i][int(root)] += header_size + math.ceil (int(count) * size / len(matrix) / flit_size)
        matrix[int(root)][i] += ack_size
    # MPI_Scatter
    elif ('MPI_Scatter entering' or 'MPI_Scatterv entering') in lineStr:
      countLine = lines[i+2].decode('utf-8')
      count = re.search('sendcount=(\d+)', countLine).group(1)
      datatypeLine = lines[i+3].decode('utf-8')
      datatype = re.search('sendtype=(\d+)', datatypeLine).group(1)
      size = MPI_Data_Type_to_size(int(datatype))
      rootLine = lines[i+6].decode('utf-8')
      root = re.search('root=(\d+)', countLine).group(1)
      for i in range(len(matrix)):
        matrix[int(root)][i] += header_size + math.ceil (int(count) * size / len(matrix) / flit_size)
        matrix[i][int(root)] += ack_size
    # MPI_Reduce
    elif ('MPI_Reduce entering') in lineStr:
      countLine = lines[i+1].decode('utf-8')
      count = re.search('count=(\d+)', countLine).group(1)
      datatypeLine = lines[i+2].decode('utf-8')
      datatype = re.search('datatype=(\d+)', datatypeLine).group(1)
      size = MPI_Data_Type_to_size(int(datatype))
      rootLine = lines[i+4].decode('utf-8')
      root = re.search('root=(\d+)', countLine).group(1)
      for i in range(len(matrix)):
        matrix[i][int(root)] += header_size + math.ceil (int(count) * size / len(matrix) / flit_size)
        matrix[int(root)][i] += ack_size
    # MPI_Allgather
    elif ('MPI_Allgather entering' or 'MPI_Allgatherv entering') in lineStr:
      countLine = lines[i+1].decode('utf-8')
      count = re.search('sendcount=(\d+)', countLine).group(1)
      datatypeLine = lines[i+2].decode('utf-8')
      datatype = re.search('sendtype=(\d+)', datatypeLine).group(1)
      size = MPI_Data_Type_to_size(int(datatype))
      for i in range(len(matrix)):
        for j in range(len(matrix[0])):
          matrix[i][j] += header_size + math.ceil (int(count) * size / len(matrix) / flit_size)
          matrix[j][i] += ack_size
    # MPI_Allreduce
    elif ('MPI_Allreduce entering') in lineStr:
      countLine = lines[i+1].decode('utf-8')
      count = re.search('count=(\d+)', countLine).group(1)
      datatypeLine = lines[i+2].decode('utf-8')
      datatype = re.search('datatype=(\d+)', datatypeLine).group(1)
      size = MPI_Data_Type_to_size(int(datatype))
      for i in range(len(matrix)):
        for j in range(len(matrix[0])):
          matrix[i][j] += header_size + math.ceil (int(count) * size / len(matrix) / flit_size)
          matrix[j][i] += ack_size
    # MPI_Reduce_Scatter
    elif ('MPI_Reduce_scatter entering') in lineStr:
      countLine = lines[i+2].decode('utf-8')
      count = re.search('recvcounts.(\d+)', countLine).group(1)
      datatypeLine = lines[i+3].decode('utf-8')
      datatype = re.search('datatype=(\d+)', datatypeLine).group(1)
      size = MPI_Data_Type_to_size(int(datatype))
      for i in range(len(matrix)):
        for j in range(len(matrix[0])):
          matrix[i][j] += header_size + math.ceil (int(count) * size / len(matrix) / flit_size)
          matrix[j][i] += ack_size
    # MPI_Put
    elif ('MPI_Put entering') in lineStr:
      countLine = lines[i+1].decode('utf-8')
      count = re.search('origincount=(\d+)', countLine).group(1)
      datatypeLine = lines[i+2].decode('utf-8')
      datatype = re.search('origintype=(\d+)', datatypeLine).group(1)
      size = MPI_Data_Type_to_size(int(datatype))
      destLine = lines[i+3].decode('utf-8')
      dest = int(re.search('targetrank=(\d+)', destLine).group(1))
      matrix[src][dest] += header_size + math.ceil (int(count) * size / flit_size)
      matrix[dest][src] += ack_size
    # MPI_Get
    elif ('MPI_Get entering') in lineStr:
      countLine = lines[i+1].decode('utf-8')
      count = re.search('origincount=(\d+)', countLine).group(1)
      datatypeLine = lines[i+2].decode('utf-8')
      datatype = re.search('origintype=(\d+)', datatypeLine).group(1)
      size = MPI_Data_Type_to_size(int(datatype))
      destLine = lines[i+3].decode('utf-8')
      dest = int(re.search('targetrank=(\d+)', destLine).group(1))
      matrix[dest][src] += header_size + math.ceil (int(count) * size / flit_size)
      matrix[src][dest] += ack_size

def main(args):
  npes = int(args.npes)
  matrix = [[0 for x in range(npes)] for y in range(npes)]
  global flit_size
  flit_size = int(args.flitsize)

  # open totol of n DUMPI logs
  for i in range(npes):
    # find all log names
    s = list(args.DUMPI)
    num = list(str(i).zfill(4))
    for j in range(4):
      s[len(s)-8+j] = num[j]
    fileName = ''.join(s)
    opener = gzip.open if fileName.endswith('.gz') else open 
    with opener(fileName, 'rb') as fd:
      analyze_DUMPI(fd, matrix)

  # compute traffic matrix
  traffic_matrix = copy.deepcopy(matrix)
  injection_rate = [0 for x in range(npes)]
  max_sum = 0
  for i in range(npes):
    sum = 0
    for j in range(len(matrix[0])):
      sum += matrix[i][j]
    injection_rate[i] = sum
    if (sum > max_sum):
      max_sum = sum
    for j in range(len(traffic_matrix[0])):
      traffic_matrix[i][j] /= sum
      traffic_matrix[i][j] = '%.4f'%(traffic_matrix[i][j])

  # compute injection rate
  for i in range(npes):
    injection_rate[i] /= max_sum
    injection_rate[i] = '%.4f'%(injection_rate[i])
  
  totaltraffic = open('totaltraffic.csv', 'w')
  writer = csv.writer(totaltraffic, delimiter=',', quotechar='"')
  writer.writerows(matrix)
  totaltraffic.close()

  trafficmatrix = open('trafficmatrix.csv', 'w')
  writer = csv.writer(trafficmatrix, delimiter=',', quotechar='"')
  writer.writerows(traffic_matrix)
  trafficmatrix.close()

  outfile = open('relativeinjection.csv', 'w')
  writer = csv.writer(outfile, delimiter=',', quotechar='"')
  for element in injection_rate:
    writer.writerow([element])
  outfile.close()


if __name__ == '__main__':
  ap = argparse.ArgumentParser()
  ap.add_argument('DUMPI',
                  help='any one of the DUMPI log file series (format: ****_000*.txt)')
  ap.add_argument('npes',
                  help='total processing elements (totol number of log files)')
  ap.add_argument('flitsize',
                  help='flit size (in bytes)')
  args = ap.parse_args()
  main(args)

