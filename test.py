import csv
import numpy as np
import pandas as pd
import Bio

from Bio.Seq import Seq

df = pd.read_csv("ABOgene.csv")

colDNA = df["Reference"]

sequenceInStr = ""
colLength = len(colDNA.index)

for i in np.arange(0, colLength, 1):
    if len(colDNA.array[i]) == 0:
        sequenceInStr = sequenceInStr + colDNA.array[i]
    else:
        sequenceInStr = sequenceInStr + colDNA.array[i][0]

dnaseq = Seq(sequenceInStr)
mRNA = dnaseq.reverse_complement().transcribe()
proteinSeq = mRNA.translate()

print(proteinSeq)

print(len(proteinSeq) * 3)
