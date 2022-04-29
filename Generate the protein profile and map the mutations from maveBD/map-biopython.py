import sys
from Bio import Entrez
from Bio import SeqIO

Entrez.email = "Your.Name.Here@example.org"

def get_mrna (prot):
	handle = Entrez.elink(dbfrom="protein", id=prot, linkname="protein_nuccore_mrna")
	record = Entrez.read(handle)
	handle.close()
	linked = [link["Id"] for link in record[0]["LinkSetDb"][0]["Link"]]
	return linked


def get_sequence (dbname,pid):
	seq=None
	handle = Entrez.efetch(db=dbname, id=pid, rettype="gbwithparts", retmode="text")
	record = SeqIO.read(handle, "genbank")
	for i in record.features:
		if i.type=='CDS':
			St = i.location.start
			End = i.location.end
			Strand = i.location.strand
			#print (St, End, Strand)
	handle.close()
	if record:
		seq=(record.id,record.seq[St:End])
	return seq




if __name__ == '__main__':
	prot=name=sys.argv[1]
	#code=prot.split('.')[0]
	mrna=get_mrna(prot)	
	#mrna=get_mrna(code)
	if len(mrna)>0:
		rna_seq=get_sequence('nucleotide',mrna[0])
		if rna_seq:
			print (">"+prot+'\n'+rna_seq[1])
			#print (">"+rna_seq[0]+'\n'+rna_seq[1])
			#prot_seq=get_sequence('protein',prot)
			#if prot_seq:
				#print (">"+prot_seq[0])
				#print (''.join([' '+aa+' ' for aa in prot_seq[1]]))
