from Bio import SeqIO
DnaB = open("/data/scratch/janani/molevolvr_out/DnaBcm_full/DnaBcm.fa", "a")
DciA = open("/data/scratch/janani/molevolvr_out/DciAcm_full/DciAcm.fa", "a")
with open("/data/scratch/janani/molevolvr_out/Dciacb_full/Dciacb.fa", "r") as in_file:
    for seq in SeqIO.parse(in_file, "fasta"):
        print(seq.id)
        if "helicase" in seq.description:
            SeqIO.write(seq, DnaB, "fasta")
        else:
            SeqIO.write(seq,DciA, "fasta" )
DnaB.close()
DciA.close()
