# -*- coding: utf-8 -*-

import numpy as np
import copy as cp

m = 9 
n = int(np.power(2,m))
T = np.array([[1, 0],
              [1, 1]]
              , dtype='int')


def recursivelyCalc_Iw(ep):
  # BEC の相互情報量を計算

  Iw = np.zeros((m+1,n))
  Iw[0][0] = 1 - ep
  for i in range(1,m+1):
    tmp = int(np.power(2,i))
    for j in range(tmp):
      if j%2 == 0:
        Iw[i][j] = np.power(Iw[i-1][j/2],2)
      else:
        Iw[i][j] = 2*Iw[i-1][int(j/2)] -  np.power(Iw[i-1][int(j/2)],2)

  return Iw[m]



def Bitrev(k):
  # 置換行列 A_k2 の作成

  k_2 = int(np.power(2, k))
  A_k2 = np.zeros((k_2, k_2), dtype='int')
  for i in range(k_2):
      if i % 2 == 0:
        tmp = int(i / 2)
        A_k2[i][tmp] = 1
      else:
        tmp = int(i / 2 + k_2 / 2)
        A_k2[i][tmp] = 1

  return A_k2


def Dot_XOR(A, B):  # ABの積(XOR)
  out = np.dot(A, B)
  for i in range(n):
    if out[i] % 2 == 0:
      out[i] = 0
    else:
      out[i] = 1
  return out



def MakeGenMatrix():
  Tn = cp.deepcopy(T)
  Bn = np.eye(2, dtype='int')
  for i in range(m - 1):
    Tn = np.kron(T, Tn)
    Bn = np.dot(Bitrev(i + 2), np.kron(np.eye(2), Bn))
  Gen = np.dot(Bn, Tn)  # 生成行列

  return Gen



def main():
  Gen = MakeGenMatrix()
  Gen = Gen.astype(np.int64)
  nameGen = 'GenerateMatrix_' + str(n) + '.txt'
  np.savetxt(nameGen, Gen, fmt='%d')


if __name__ == '__main__':
    main()
