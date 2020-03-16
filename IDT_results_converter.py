#!/usr/loca/env python3

import sys
import pandas as pd

idt_file = sys.argv[1]
fasta_out = sys.argv[2]
prefix = sys.argv[3]

df = pd.read_excel(idt_file)

type_dict = {
    'Forward Primer': '-F',
    'Probe': '-P',
    'Reverse Primer': '-R'
}

with open(fasta_out, 'w') as f:
    i = 0
    assay_list = list()
    for index, row in df.iterrows():
        typ = row['Type']
        seq = row['Sequence']
        size = row['Amplicon']

        assay_list.append((typ, seq, size))
        if len(assay_list) == 4:
            j = i
            if prefix:
                i = '{}_{}'.format(prefix, i)
            f.write('>{}_{}bp{}\n{}\n'.format(i, int(assay_list[3][2]), type_dict[assay_list[0][0]], assay_list[0][1]))
            f.write('>{}_{}bp{}\n{}\n'.format(i, int(assay_list[3][2]), type_dict[assay_list[1][0]], assay_list[1][1]))
            f.write('>{}_{}bp{}\n{}\n'.format(i, int(assay_list[3][2]), type_dict[assay_list[2][0]], assay_list[2][1]))
            f.write('\n')
            assay_list = list()
            i = j
            i += 1
