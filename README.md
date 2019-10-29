# CASCL
CRC-aided successive cancellation list decoder of polar codes

- Construction of Polar Codes
I use bhattacharyya analisis to select frozen bit sets.
PolarMatrix.py gives the information of generate matrix G (for example GenerateMatrix_2048.txt). 

- Encode  
x = u G
G is sparse matrix(I use COO), so it takes O(N log N) time complexity.

- Decode
CRC-16 (you can change the length)  
BPSK & AWGN channel
