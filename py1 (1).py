from Bio import SeqIO
for record in SeqIO.parse("hp_gene.fasta", "fasta"):
    print('ID:',record.id)
    print('Description:',record.description)
    print('Sequence:',record.seq)   
    print ('Length:',len(record.seq))   
    print("Annotations:",record.annotations )
    print(len(record.features))
    print(record.features)


#step 2: Sequence quality analysis
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

record = SeqIO.read("hp_gene.fasta", "fasta")

sequence_length = len(record.seq)
gc_content = gc_fraction(record.seq) * 100

print("Sequence ID:", record.id)
print("Sequence Length:", sequence_length)
print("GC Content (%):", gc_content)


#step 3: sequence filtering and validation
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
if sequence_length >= 300 and 30 <= gc_content <= 70:
    print("Sequence passed quality filtering")
    SeqIO.write(record, "filtered_hp_gene.fasta", "fasta")
else:
    print("Sequence failed quality filtering")

# #step 4: sequence alignment using BLAST
from Bio.Blast import NCBIWWW
result_handle = NCBIWWW.qblast(
    program="blastn",
    database="nt",  # nucleotide collection database
    sequence=record.seq
    )
with open("blastn_result.xml", "w") as b:  
    b.write(result_handle.read())  
print("BLAST search completed successfully")  

from Bio.Blast import NCBIXML
with open("blastn_result.xml") as r:
    blast_record = NCBIXML.read(r)

print('Number of alignments:', len(blast_record.alignments))
best_alignment = blast_record.alignments[:10]    

for alignment in best_alignment:
    print('Alignment title:', alignment.title)
    print('Length:', alignment.length)
    for hsp in alignment.hsps:
        print('E-value:', hsp.expect)
        print("Score:", hsp.score)
        print("Query sequence:", hsp.query)
        print("Alignment match:", hsp.match)
        print("Matched sequence:", hsp.sbjct)
        print("Query range:", hsp.query_start, "-", hsp.query_end)
        print("Subject range:", hsp.sbjct_start, "-", hsp.sbjct_end)
    print("--"*50) 

print('Evolutionary hints found in:')

count = 1
for alignment in best_alignment:
    if any(hsp.expect < 0.01 for hsp in alignment.hsps):
        print('Alignment title',count, ':', alignment.title)
        count += 1
           
