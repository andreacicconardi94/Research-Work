#!/home/blasfemus/anaconda3/bin/python

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
from Bio import SeqIO
import sys

Entrez.email = "Your.Name.Here@example.org"

#THIS CODE WILL DOWNLOAD AN ALIGNMENT FILE FROM BLASTP
#THEN IT WILL PARSE IT AND PRODUCE THE ALIGNMENT FILE
#AS STABLE_ALIGNMENT.FASTA

def Blastp(query):
	fastafile = open(query).readlines()[1:]
	prot_seq = ''
	for line in fastafile:
		line.strip()
		prot_seq += line

	'''
	result_handle = NCBIWWW.qblast("blastp", "refseq_protein", prot_seq, expect=0.05, hitlist_size=500)
	with open("my_blast.xml", "w") as out_handle:
		out_handle.write(result_handle.read())
		blast_records = NCBIXML.parse(result_handle)

	result_handle.close()
	'''

	#SET THE PATH OF my_blast.xml, IT SHOULD BE THE SAME OF THE CODE FOLDER
	out_handle = open(r'/home/blasfemus/Desktop/Nuova_cartella/Lavoro/Synonymous/Working/my_blast.xml')
	blast_record = NCBIXML.parse(out_handle)

	id_list = []
	id_dict = {}
	for record in blast_record:
		for align in record.alignments:
			for hsp in align.hsps:
				id_dict[align.hit_id[4:-1]] = str(hsp.sbjct_start)+'-'+str(hsp.sbjct_end)
				id_list.append(align.hit_id[4:-1])

	#THIS WILL PRODUCE A FILE CONTAINING THE HEADERS OF PROTEIN
	id_list = list(set(id_list))
	header_file = open("list_of_header.txt", 'w')
	for id_prot in id_list:
		header_file.write(id_prot+'\n')

	header_file.close()


	stable_align = open("alignment.fasta", 'w')
	
	for identifier in id_dict.keys():
		handle = Entrez.efetch(db="protein", id=identifier, rettype="fasta", retmode="txt")
		record = SeqIO.read(handle, "fasta")
		stable_align.write('>'+record.id+':'+id_dict[record.id]+'\n'+str(record.seq)+'\n')

	stable_align.close()



#THE INPUT PROTEIN SEQUENCE CAN BE INSERED WITH ARGV FUNCTION
#OR WRITING THE WHOLE PATH OF THE FASTA FILE

if __name__=='__main__':
	Protein = r'/home/blasfemus/Desktop/Nuova_cartella/Lavoro/Synonymous/Working/protein_sequence_01-b-2.fasta'
	#Protein = sys.argv[1]
	Blastp(Protein)
