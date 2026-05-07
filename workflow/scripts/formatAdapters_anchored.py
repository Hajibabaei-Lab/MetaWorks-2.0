import argparse

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Create anchored linked adapter FASTA from a CSV file."
    )
    parser.add_argument("filename", help="CSV file with adapter definitions")
    args = parser.parse_args(argv)

    file_out = 'adapters.fasta'

    df = pd.read_csv(args.filename)

    df['SampleIDamplicon'] = df['SampleID'].str.cat(df['Amplicon'], sep="_")
    df['SampleIDamplicon'] = df['SampleIDamplicon'].astype(str) + ";"

    df['Forward'] = df['Forward'].str.replace('I', 'N')
    df['Reverse'] = df['Reverse'].str.replace('I', 'N')

    i = 0
    rc = [None] * len(df['Reverse'])
    for x in df['Reverse']:
        seq = Seq(x)
        rc[i] = str(seq.reverse_complement())
        i = i + 1

    df['ReverseRC'] = rc

    df['LinkedAdapters'] = df['Forward'].str.cat(df['ReverseRC'], sep="...")
    df['LinkedAdapters'] = df['LinkedAdapters'].map('^{}$'.format)

    record_list = []
    with open(file_out, 'w') as f_out:
        for index, row in df.iterrows():
            record = SeqRecord(
                Seq(df['LinkedAdapters'].iloc[index]),
                description="",
                id=df['SampleIDamplicon'].iloc[index],
            )
            record_list.append(record)
        SeqIO.write(record_list, f_out, 'fasta')


if __name__ == "__main__":
    main()
