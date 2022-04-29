from Bio.Data import CodonTable
import pandas as pd


def get_codon():
	c=CodonTable.standard_dna_table.forward_table
	d={}
	for i in c.keys():
		d[c[i]]=d.get(c[i],[])
		d[c[i]].append(i)

	return d


def get_degeneri(profile, AA, COD, dizionario):
	Degeneri = dizionario[AA]
	fCOD = profile.get(COD, 0)
	Frequencies = {}
	for i in Degeneri:
		if i != COD:
			f = profile.get(i, 0)
			if f > 0:
				Frequencies[i] = f

	return fCOD, Frequencies




def get_non_synonimous(profile, Frequencies, COD):
	NS_dic = {}

	for j in Frequencies:
		profile.pop(j, None)

	for k in profile:
		if k != COD and k != 'name':
			NS_dic[k] = profile[k]


	return  NS_dic




def parse_prot(filepro, fileseq):
	protseq = filepro.read().splitlines()
	Startlist = []
	Endlist = []
	Sequences = []
	HEADER = []
	for i in protseq:
		if i[0]=='>':
			test = i.split(':')
			positions = test[1].split()
			end = int(positions[0].split('-')[1])
			start = int(positions[0].split('-')[0])
			Startlist.append(start)
			Endlist.append(end)
			HEADER.append(i)
		else:
			Sequences.append(i)

	rnaseq = fileseq.read().splitlines()
	RNA = []
	for a in rnaseq:
		if a[0]=='>':
			continue
		else:
			RNA.append(a)

	nali=len(Sequences[0])
	count = {}
	Codons = {}
	c = 0
	for b in range(0, nali):
		if Sequences[0][b] != '-':
			count[b] = c
			Codons[c] = (Sequences[0][b], RNA[0][(Startlist[0]-1+c)*3:(Startlist[0]-1+c)*3+3])
			c += 1
		else: continue
	#print (Codons)


	Profile = {}
	for x in range(nali):
		Profile[x] = {}

	for j in range(0,len(Sequences)):
		#print (Startlist[j])
		c = Startlist[j]-1
		Profile = update_profile(Profile, Sequences[j], RNA[j], c)

		'''
		for k in range(nali):
			if Sequences[j][k] == '-':
				continue
			else:
				codon = RNA[j][c*3:c*3+3]
				if Profile.get(k, 0) != 0:
					Profile[k][codon] = Profile[k].get(codon, 0) + 1
				c += 1
		'''

	Profile1 = {}
	for i in count.keys():
		Profile1[count[i]] = Profile[i]
		Profile1[count[i]]['name'] = Codons[count[i]]

	return Profile1


def update_profile (profile, protein, rna, startlist):
	c = startlist
	nali = len(protein)
	for k in range(nali):
		if protein[k] == '-':
			continue
		else:
			codon = rna[c*3:c*3+3]
			if profile.get(k, 0) != 0:
				profile[k][codon] = profile[k].get(codon, 0) + 1
			c += 1

	return profile




if __name__=='__main__':
	Protein = open(r'/home/blasfemus/Desktop/Nuova_cartella/Lavoro/Synonymous/alignment_protein_01-b-2.fasta')
	DNA = open(r'/home/blasfemus/Desktop/Nuova_cartella/Lavoro/Synonymous/stable_risultati_prova_01-b-2.txt')
	dizionario = get_codon()
	Profile = parse_prot(Protein, DNA)
	n = list(Profile.keys())
	n.sort()

	list_k = []
	list_aa = []
	list_cod = []
	list_fcod = []
	list_i = []
	list_freqdeg = []
	list_ndeg = []
	list_ntotal = []


	ns_list_k = []
	ns_list_aa = []
	ns_list_cod = []
	ns_list_fcod = []
	ns_list_j = []
	ns_list_freqcodonenonsinonimo = []
	ns_list_n_ns = []
	ns_list_ntotal = []


	for k in n:
		aa = Profile[k]['name'][0]
		cod = Profile[k]['name'][1]
		df = Profile[k]
		del(df['name'])
		ntotal = sum(list(df.values()))


		fcod, freqdegeneri = get_degeneri(Profile[k], aa, cod, dizionario)
		NS_dictionary = get_non_synonimous(Profile[k], freqdegeneri, cod)


		ndeg = sum(list(freqdegeneri.values()))
		n_ns = sum(list(NS_dictionary.values()))


		for i in freqdegeneri.keys():
			#print (k, aa, cod, fcod, i, freqdegeneri[i], ndeg, ntotal)
			list_k.append(k)
			list_aa.append(aa)
			list_cod.append(cod)
			list_fcod.append(fcod)
			list_i.append(i)
			list_freqdeg.append(freqdegeneri[i])
			list_ndeg.append(ndeg)
			list_ntotal.append(ntotal)



		for j in NS_dictionary.keys():
			ns_list_k.append(k)
			ns_list_aa.append(aa)
			ns_list_cod.append(cod)
			ns_list_fcod.append(fcod)
			ns_list_j.append(j)
			ns_list_freqcodonenonsinonimo.append(NS_dictionary[j])
			ns_list_n_ns.append(n_ns)
			ns_list_ntotal.append(ntotal)



	table = {'Posizione' : list_k,
			'Residuo' : list_aa,
			'Codone_originale' : list_cod,
			'Frequenza_originale' : list_fcod,
			'Codone_degenere' : list_i,
			'Frequenza_degenere' : list_freqdeg,
			'Numero_degeneri' : list_ndeg,
			'Numero_totali' : list_ntotal
			}


	NS_table = {'Posizione' : ns_list_k,
			'Residuo' : ns_list_aa,
			'Codone_originale' : ns_list_cod,
			'Frequenza_degenere' : ns_list_fcod,
			'Codone_non_sinonimo' : ns_list_j,
			'Frequenza_codone_non_sinonimo' : ns_list_freqcodonenonsinonimo,
			'Numero_totale_non_sinonimo' : ns_list_n_ns,
			'Numero_totale_codoni' : ns_list_ntotal
			}

	df = pd.DataFrame(table)

	ns_df = pd.DataFrame(NS_table)

	df.to_csv(r'/home/blasfemus/Desktop/Nuova_cartella/Lavoro/Synonymous/Working/Tabella_degeneri_01-b-2.csv')
	ns_df.to_csv(r'/home/blasfemus/Desktop/Nuova_cartella/Lavoro/Synonymous/Working/Tabella_non_sinonimi_01-b-2.csv')


		#print (k, aa, cod, fcod, freqdegeneri, sum(list(df.values())), sum(list(freqdegeneri.values())))
