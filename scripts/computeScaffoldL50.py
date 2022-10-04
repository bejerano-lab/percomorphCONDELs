#!/cluster/u/hchen17/miniconda/bin/python

import sys
import pandas as pd

chromSizes = pd.read_csv(sys.argv[1], sep='\t', header=None)

# print(chromSizes.columns[1])

chromSizes = chromSizes.sort_values(chromSizes.columns[1], ascending=False).reset_index(drop=True)

# stream = os.popen('sort -k2,2nr '+chromSizes)
# output = stream.read()

assemblyLength = chromSizes[1].sum()
halfLength = int(assemblyLength/2)

# print(chromSizes.head())
# print(assemblyLength)

L50 = 0
subtotal = 0

while subtotal < halfLength:
	L50 += 1
	subtotal += chromSizes[1][L50-1]
	# print("subtotal: " + str(subtotal))

print(L50)