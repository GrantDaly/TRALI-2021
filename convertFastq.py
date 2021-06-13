from Bio import SeqIO
from Bio.Seq import Seq

file_out = "somatic-masked.fa"
with open(file_out, 'w') as f_out:
    for seq_record in SeqIO.parse("../ensembl-98/Homo_sapiens.GRCh38.dna.primary_assembly.ens98.chr.fa", "fasta"):
        print(seq_record.id)
        tempLength = len(seq_record)
        print(tempLength)
        if(seq_record.id != "chrM"):
            tempSeqObject = Seq("N" * tempLength)
            seq_record.seq = tempSeqObject
        SeqIO.write(seq_record, f_out, 'fasta')
