# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 13:32:53 2021

@author: Andrea
"""

import pandas as pd
from scipy.stats import pearsonr
import math
from Bio.SeqUtils import seq3
import numpy as np
import matplotlib.pyplot as plt




#FIRST FUNCTION TO EXTRACT INFORMATION FROM MAVEDB
def extract_mutations(Mutations_file):
    Muts_file = Mutations_file.readlines()
    Positions = []
    PositionsAA = []
    TempMutCod = []
    AA = []
    Scores = []
    
    NS_Positions = []
    NS_PositionsAA = []
    NS_TempMutCod = []
    NS_AA = []
    NS_Scores = []
    

    for mutations in Muts_file:
        mutations = mutations.strip()
        mutations = mutations.split(' ')
        try:
            posAA = int(mutations[1][5:-1])-1
            pos = int(mutations[0][0:-3])
            changed_nucleotide = mutations[0][-1]
            Positions.append(pos)
            PositionsAA.append(posAA)
            TempMutCod.append(changed_nucleotide)
            Residue = mutations[1][2:5]
            AA.append(Residue)
            score = float(mutations[2])
            Scores.append(score)

        except:
            NS_posAA = int(mutations[1][5:-3])-1
            pos = int(mutations[0][0:-3])
            changed_nucleotide = mutations[0][-1]
            NS_Positions.append(pos)
            NS_PositionsAA.append(NS_posAA)
            NS_TempMutCod.append(changed_nucleotide)
            Residue = mutations[1][2:5]
            NS_AA.append(Residue)
            score = float(mutations[2])
            NS_Scores.append(score)



    return (Positions, TempMutCod, AA, Scores, PositionsAA, NS_Positions, NS_TempMutCod, NS_AA, NS_Scores, NS_PositionsAA)



#SECOND FUNCTION TO EXTRACT SEQUENCE
def extract_sequence(Sequence):
    next(Sequence)
    Seq = Sequence.readlines()
    sequence_dna = ''
    sequence_rna = ''
    for line in Seq:
        line = line.strip()
        sequence_dna += line

    for i in sequence_dna:
        if i == "T":
            sequence_rna += "U"
        else:
            sequence_rna += i

    return (sequence_rna, sequence_dna)




#THIRD FUNCTION TO EXTRACT THE CODONS AND THE MUTANTS
def extract_statistics(sequence_rna, sequence_dna, TempMutCod, Positions):
    WT_codons = []
    Mut_codons = []
    for i in range(len(Positions)):
        if Positions[i]%3 == 0:
            codon = sequence_dna[Positions[i]-3:Positions[i]]
            WT_codons.append(codon)
            Mutant = codon[0]+codon[1]+TempMutCod[i]
            Mut_codons.append(Mutant)
        elif Positions[i]%3 == 1:
            codon = sequence_dna[Positions[i]-1:Positions[i]+2]
            WT_codons.append(codon)
            Mutant = TempMutCod[i]+codon[1]+codon[2]
            Mut_codons.append(Mutant)
        elif Positions[i]%3 == 2:
            codon = sequence_dna[Positions[i]-2:Positions[i]+1]
            WT_codons.append(codon)
            Mutant = codon[0]+TempMutCod[i]+codon[2]
            Mut_codons.append(Mutant)
        else:
            print ("ERROR")
            print (Positions[i]%3)

    '''
    for x in range (0,len(WT_codons)):
        Mutant = WT_codons[x][0]+WT_codons[x][1]+TempMutCod[x]
        Mut_codons.append(Mutant)

    Mut_codons = []
    '''
    
    Freq_mut = []
    Freq_wt = []

    Freq_organism_wt = []
    Freq_organism_mut = []

    '''
    for x in range (0,len(WT_codons)):
        if TempMutCod[x] == 'T':
            Mutant = WT_codons[x][0]+WT_codons[x][1]+'U'
        else:
            Mutant = WT_codons[x][0]+WT_codons[x][1]+TempMutCod[x]

        Mut_codons.append(Mutant)
    '''
    
    
    Homo_sapiens = create_codon_table()
    Saccharomyces_cerevisiae = create_codon_table()
    
    for x in range (0, len(WT_codons)):
        Freq_wt.append(Homo_sapiens[WT_codons[x]])
        Freq_mut.append(Homo_sapiens[Mut_codons[x]])
        Freq_organism_wt.append(Saccharomyces_cerevisiae[WT_codons[x]])
        Freq_organism_mut.append(Saccharomyces_cerevisiae[Mut_codons[x]])
    




    return (Mut_codons, Freq_wt, Freq_mut, Freq_organism_wt, Freq_organism_mut, WT_codons)




def Table_and_Corr(Positions, AA, WT_codons, Mut_codons, Freq_wt, Freq_mut, Freq_organism_wt, Freq_organism_mut, Scores, Degeneri,
                   NS_Positions, NS_AA, NS_WT_codons, NS_Mut_codons, NS_Scores, Non_sinonimi):
					   
    Table = { 'Protein' : ['UBE2I' for i in range(0, len(Positions))],
             'Position': Positions,
             'AA' : AA,
             'WT_codon': WT_codons,
             'Mut_codon': Mut_codons,
             #'FreqHumWt': Freq_wt,
             #'FreqHumMut': Freq_mut,
             #'FreqOrgWT' : Freq_organism_wt,
             #'FreqOrgMut': Freq_organism_mut,
             'Score': Scores
             }

    df = pd.DataFrame(Table)
    tri = df.sort_values(by='Position')
    Degeneri.drop(Degeneri.columns[0], 1)

    Frames = [tri, Degeneri]
    Final_df = pd.concat(Frames, axis=1)
    #Final_df1 = Final_df.sort_values(by=['Position','Posizioni'])

    Final_df.to_csv(r'/home/blasfemus/Desktop/Nuova_cartella/Lavoro/Synonymous/Final_df_01-b-2.csv')


    Degeneri_list = dataframe_to_list(Degeneri)
    
    Posizione = Degeneri_list[0]
    Residuo = Degeneri_list[1]
    Codone_originale = Degeneri_list[2]
    Frequenza_originale = Degeneri_list[3]
    Codone_degenere = Degeneri_list[4]
    Frequenza_degenere = Degeneri_list[5]
    Numero_degeneri = Degeneri_list[6]
    Numero_totali = Degeneri_list[7]

    Residuo = convert_residues3to1(Residuo)


    MaveDB_list = dataframe_to_list(tri)

    Positions = MaveDB_list[1]
    AA = MaveDB_list[2]
    WT_codons = MaveDB_list[3]
    Mut_codons = MaveDB_list[4]
    Scores = MaveDB_list[5]




    Non_sinonimi_list = dataframe_to_list(Non_sinonimi)
    
    
    NS_Posizione = Non_sinonimi_list[0]
    NS_Residuo = Non_sinonimi_list[1]
    NS_Codone_originale = Non_sinonimi_list[2]
    NS_Frequenza_originale = Non_sinonimi_list[3]
    Codone_non_sinonimo = Non_sinonimi_list[4]
    Frequenza_non_sinonimo = Non_sinonimi_list[5]
    Numero_non_sinonimi = Non_sinonimi_list[6]
    NS_Numero_totali = Non_sinonimi_list[7]

    NS_Residuo = convert_residues3to1(NS_Residuo)

    Information_dict = {}
    NS_Information_dict = {}



    for b in range(len(NS_Positions)):
        if NS_Information_dict.get((NS_Positions[b], NS_AA[b]), 0) == 0:
            NS_Information_dict[(NS_Positions[b], NS_AA[b])] = {(NS_WT_codons[b], NS_Mut_codons[b]): {'Score': NS_Scores[b]}}
        else:
            new_key = {(NS_WT_codons[b], NS_Mut_codons[b]): {'Score': NS_Scores[b]}}
            NS_Information_dict[(NS_Positions[b], NS_AA[b])].update(new_key)



    for c in range(len(Positions)):
        if Information_dict.get((Positions[c], AA[c]), 0) == 0:
            Information_dict[(Positions[c], AA[c])] = {(WT_codons[c], Mut_codons[c]): {'Score': Scores[c]}}
        else:
            new_key = {(WT_codons[c], Mut_codons[c]): {'Score': Scores[c]}}
            Information_dict[(Positions[c], AA[c])].update(new_key)



    for d in range(len(Posizione)):
        if Information_dict.get((Posizione[d], Residuo[d]), 0) != 0:
            G1 = list(Information_dict[(Posizione[d], Residuo[d])].keys())[0][0]
            G2 = list(Information_dict[(Posizione[d], Residuo[d])].keys())[0][1]
                
            if Information_dict[(Posizione[d], Residuo[d])].get((Codone_originale[d], Codone_degenere[d]), 0) != 0:
                Information_dict[(Posizione[d], Residuo[d])][(Codone_originale[d], Codone_degenere[d])]['Frequenza originale'] = Frequenza_originale[d]
                Information_dict[(Posizione[d], Residuo[d])][(Codone_originale[d], Codone_degenere[d])]['Frequenza degenere'] = Frequenza_degenere[d]
                Information_dict[(Posizione[d], Residuo[d])][(Codone_originale[d], Codone_degenere[d])]['Numero degeneri'] = Numero_degeneri[d]
                Information_dict[(Posizione[d], Residuo[d])][(Codone_originale[d], Codone_degenere[d])]['Numero totali'] = Numero_totali[d]
                
            elif Information_dict[(Posizione[d], Residuo[d])][(G1, G2)].get('Frequenza degenere', 0) == 0:
                
                if list(Information_dict[(Posizione[d], Residuo[d])].keys())[0][0] == Codone_originale[d]:
                    for m in range(len(list(Information_dict[(Posizione[d], Residuo[d])].keys()))):
                        Mutante = list(Information_dict[(Posizione[d], Residuo[d])].keys())[m][1]
                        Information_dict[(Posizione[d], Residuo[d])][(Codone_originale[d], Mutante)]['Frequenza originale'] = Frequenza_originale[d]
                        Information_dict[(Posizione[d], Residuo[d])][(Codone_originale[d], Mutante)]['Frequenza degenere'] = 0
                        Information_dict[(Posizione[d], Residuo[d])][(Codone_originale[d], Mutante)]['Numero degeneri'] = Numero_degeneri[d]
                        Information_dict[(Posizione[d], Residuo[d])][(Codone_originale[d], Mutante)]['Numero totali'] = Numero_totali[d]



    for e in range(len(NS_Posizione)):
        if NS_Information_dict.get((NS_Posizione[e], NS_Residuo[e]), 0) != 0:
            G1 = list(NS_Information_dict[(NS_Posizione[e], NS_Residuo[e])].keys())[0][0]
            G2 = list(NS_Information_dict[(NS_Posizione[e], NS_Residuo[e])].keys())[0][1]
                
            if NS_Information_dict[(NS_Posizione[e], NS_Residuo[e])].get((NS_Codone_originale[e], Codone_non_sinonimo[e]), 0) != 0:
                NS_Information_dict[(NS_Posizione[e], NS_Residuo[e])][(NS_Codone_originale[e], Codone_non_sinonimo[e])]['Frequenza originale'] = NS_Frequenza_originale[e]
                NS_Information_dict[(NS_Posizione[e], NS_Residuo[e])][(NS_Codone_originale[e], Codone_non_sinonimo[e])]['Frequenza non sinonimo'] = Frequenza_non_sinonimo[e]
                NS_Information_dict[(NS_Posizione[e], NS_Residuo[e])][(NS_Codone_originale[e], Codone_non_sinonimo[e])]['Numero non sinonimi'] = Numero_non_sinonimi[e]
                NS_Information_dict[(NS_Posizione[e], NS_Residuo[e])][(NS_Codone_originale[e], Codone_non_sinonimo[e])]['Numero totali'] = NS_Numero_totali[e]
                
            elif NS_Information_dict[(NS_Posizione[e], NS_Residuo[e])][(G1, G2)].get('Frequenza non sinonimo', 0) == 0:
				
                if list(NS_Information_dict[(NS_Posizione[e], NS_Residuo[e])].keys())[0][0] == NS_Codone_originale[e]:
                    for m in range(len(list(NS_Information_dict[(NS_Posizione[e], NS_Residuo[e])].keys()))):
                        Mutante = list(NS_Information_dict[(NS_Posizione[e], NS_Residuo[e])].keys())[m][1]
                        NS_Information_dict[(NS_Posizione[e], NS_Residuo[e])][(NS_Codone_originale[e], Mutante)]['Frequenza originale'] = NS_Frequenza_originale[e]
                        NS_Information_dict[(NS_Posizione[e], NS_Residuo[e])][(NS_Codone_originale[e], Mutante)]['Frequenza non sinonimo'] = 0
                        NS_Information_dict[(NS_Posizione[e], NS_Residuo[e])][(NS_Codone_originale[e], Mutante)]['Numero non sinonimi'] = Numero_non_sinonimi[e]
                        NS_Information_dict[(NS_Posizione[e], NS_Residuo[e])][(NS_Codone_originale[e], Mutante)]['Numero totali'] = NS_Numero_totali[e]





    DeltaFrqHum = []
    DeltaFrqOrganism =[]

    write_file(Information_dict, 'Risultati_dict.txt')
    write_file_NS(NS_Information_dict, 'Risultati_dict_NS.txt')
    
    Sc, Fr_wt, Fr_deg, Num_deg, Num_tot = extract_from_file(Information_dict)
    
    NS_Sc, NS_Fr_wt, NS_Fr_deg, NS_Num_deg, NS_Num_tot = extract_from_file_NS(NS_Information_dict)
    

    for i in range(0, len(Freq_wt)):
        DeltaTempHum = math.sqrt(math.pow(Freq_wt[i]-Freq_mut[i], 2))
        DeltaTempOrg = math.sqrt(math.pow(Freq_organism_wt[i]-Freq_organism_mut[i], 2))
    
        DeltaFrqHum.append(DeltaTempHum)
        DeltaFrqOrganism.append(DeltaTempOrg)


    corr1, _ = pearsonr(Freq_wt, Scores)
    corr2, _ = pearsonr(Scores, Freq_mut)
    corr3, _ = pearsonr(Scores, Freq_organism_wt)
    corr4, _ = pearsonr(Scores, Freq_organism_mut)
    corr5, _ = pearsonr(Scores, DeltaFrqHum)
    corr6, _ = pearsonr(Scores, DeltaFrqOrganism)


    print('Pearsons correlation Scores-FreqHumWT: %.3f' % corr1)
    print('Pearsons correlation Scores-FreqHumMut: %.3f' % corr2)
    print('Pearsons correlation Scores-FreqOrganismWT: %.3f' % corr3)
    print('Pearsons correlation Scores-FreqOrganismMut: %.3f' % corr4)
    print('Pearsons correlation Scores-DeltaFreqHumMut: %.3f' % corr5)
    print('Pearsons correlation Scores-DeltaFreqOrgMut: %.3f' % corr6)



    NUM_WT = np.array(Fr_wt)
    NUM_MUT = np.array(Fr_deg)
    NUM_DEG = np.array(Num_deg)
    NUM_TOT = np.array(Num_tot)

    SC_array = np.array(Sc)


    NS_NUM_WT = np.array(NS_Fr_wt)
    NS_NUM_MUT = np.array(NS_Fr_deg)
    NS_NUM_DEG = np.array(NS_Num_deg)
    NS_NUM_TOT = np.array(NS_Num_tot)
    
    NS_SC_array = np.array(NS_Sc)



    FR_WT = NUM_WT/NUM_TOT
    FR_MUT = NUM_MUT/NUM_TOT
    FR_DEG = NUM_DEG/NUM_TOT
    R1 = NUM_DEG/NUM_WT
    R2 = NUM_MUT/NUM_WT
    R3 = (NUM_TOT)

    NS_FR_WT = NS_NUM_WT/NS_NUM_TOT
    NS_FR_MUT = NS_NUM_MUT/NS_NUM_TOT
    NS_FR_DEG = NS_NUM_DEG/NS_NUM_TOT
    NS_R1 = NS_NUM_DEG/NS_NUM_WT
    NS_R2 = NS_NUM_MUT/NS_NUM_WT
    NS_R3 = (NS_NUM_TOT)

    index = []
    for x in range(len(NS_R2)):
        if NS_R2[x] == 0.0:
            index.append(x)


    NS_R2_new = np.delete(NS_R2, index)
    NS_SC_array_new = np.delete(NS_SC_array, index)


    corr1, _ = pearsonr(SC_array, FR_WT)
    corr2, _ = pearsonr(SC_array, FR_MUT)
    corr3, _ = pearsonr(SC_array, FR_DEG)
    corr4, _ = pearsonr(SC_array, R1)
    corr5, _ = pearsonr(SC_array, R2)
    #corr6, _ = pearsonr(Sc, (Num_deg/Num_tot))

    #plot_results(SC_array, NUM_WT, 'Frequenza Wyld Type')
    #plot_results(R1, SC_array, 'Numero Degeneri su Numero WT')
    #plot_results(R2, SC_array, 'Numero Mutanti su Numero WT')
    
    plot_results(NS_R1, R1, NS_SC_array, SC_array, 'Numero totale di NS su Numero totale di WT', 'Numero totale Degeneri su Numero totale di WT')
    plot_results(NS_R2_new, R2, NS_SC_array_new, SC_array, 'Numero di NS di quel tipo su Numero totale di WT', 'Numero di Degenere di quel tipo su Numero totale di WT')
    
    #plot_results(SC_array, FR_MUT, 'Frequenza Mutante')
    #plot_results(SC_array, FR_DEG, 'Frequenza Degeneri')


    print('Pearsons correlation Scores-Frequency_WT: %.3f' % corr1)
    print('Pearsons correlation Scores-Frequency_Mut: %.3f' % corr2)
    print('Pearsons correlation Scores-Number_Deg: %.3f' % corr3)
    print('Pearsons correlation Scores-Ratio1: %.3f' % corr4)
    print('Pearsons correlation Scores-Ratio2: %.3f' % corr5)
    #print('Pearsons correlation Scores-Ratio_Number: %.3f' % corr6)





def plot_results(X, Z, Y, W, X_label, Z_label):
    fig = plt.figure()
    ax1 = fig.add_subplot(111) 
    ax1.scatter(X, Y, s=10, c='b', marker="s", label=X_label)
    ax1.scatter(Z, W, s=10, c='r', marker="o", label=Z_label)
    plt.legend(loc='upper left');
    plt.show()




def write_file(dictionary, filename):
    file = open(r'/home/blasfemus/Desktop/Nuova_cartella/Lavoro/Synonymous/Working' + filename, 'w')
    Header = '#AA'+'\t'+'POS'+'\t'+'COD_WT'+'\t'+'COD_MUT'+'\t'+'MAVEDB'+'\t'+'N_WT'+'\t'+'N_DEG'+'\t'+'N_DEGS'+'\t'+'N_TOT'+'\n'
    file.write(Header)
    for key in dictionary:
        for sub_key in dictionary[key]:
            RES = str(key[1])
            POS = str(key[0])
            WT_C = str(sub_key[0])
            MUT_C = str(sub_key[1])
            SCORE = str(dictionary[key][sub_key]['Score'])
            try:
                N_WT = str(dictionary[key][sub_key]['Frequenza originale'])
                N_DEG = str(dictionary[key][sub_key]['Frequenza degenere'])
                N_DEGS = str(dictionary[key][sub_key]['Numero degeneri'])
                N_TOT = str(dictionary[key][sub_key]['Numero totali'])
            except:
                print (dictionary[key])
                print (key)
                print ('Error')
                continue

            string = RES+'\t'+POS+'\t'+WT_C+'\t'+MUT_C+'\t'+SCORE+'\t'+N_WT+'\t'+N_DEG+'\t'+N_DEGS+'\t'+N_TOT+'\n'
            file.write(string)
    
    file.close()



def write_file_NS(dictionary, filename):
    file = open(r'/home/blasfemus/Desktop/Nuova_cartella/Lavoro/Synonymous/Working' + filename, 'w')
    Header = '#AA'+'\t'+'POS'+'\t'+'COD_WT'+'\t'+'COD_NS'+'\t'+'MAVEDB'+'\t'+'N_WT'+'\t'+'N_NS'+'\t'+'N_TOT_NS'+'\t'+'N_TOT'+'\n'
    file.write(Header)
    for key in dictionary:
        for sub_key in dictionary[key]:
            RES = str(key[1])
            POS = str(key[0])
            WT_C = str(sub_key[0])
            MUT_C = str(sub_key[1])
            SCORE = str(dictionary[key][sub_key]['Score'])
            try:
                N_WT = str(dictionary[key][sub_key]['Frequenza originale'])
                N_DEG = str(dictionary[key][sub_key]['Frequenza non sinonimo'])
                N_DEGS = str(dictionary[key][sub_key]['Numero non sinonimi'])
                N_TOT = str(dictionary[key][sub_key]['Numero totali'])
            except:
                #print (dictionary[key])
                #print (key)
                #print ('Error_NS')
                continue

            string = RES+'\t'+POS+'\t'+WT_C+'\t'+MUT_C+'\t'+SCORE+'\t'+N_WT+'\t'+N_DEG+'\t'+N_DEGS+'\t'+N_TOT+'\n'
            file.write(string)
    
    file.close()


def extract_from_file(dictionary):
    Score = []
    Freq_wt = []
    Freq_deg = []
    N_degs = []
    N_tot = []
    
    for key in dictionary:
        for sub_key in dictionary[key]:
            SCORE = dictionary[key][sub_key]['Score']
            try:
                FREQ_WT = dictionary[key][sub_key]['Frequenza originale']
                FREQ_DEG = dictionary[key][sub_key]['Frequenza degenere']
                N_DEGS = dictionary[key][sub_key]['Numero degeneri']
                N_TOT = dictionary[key][sub_key]['Numero totali']
                
                Score.append(SCORE)
                Freq_wt.append(FREQ_WT)
                Freq_deg.append(FREQ_DEG)
                N_degs.append(N_DEGS)
                N_tot.append(N_TOT)
                
            except:
                continue
                      
    
    return Score, Freq_wt, Freq_deg, N_degs, N_tot
    


def extract_from_file_NS(dictionary):
    Score = []
    Freq_wt = []
    Freq_deg = []
    N_degs = []
    N_tot = []
    
    for key in dictionary:
        for sub_key in dictionary[key]:
            SCORE = dictionary[key][sub_key]['Score']
            try:
                FREQ_WT = dictionary[key][sub_key]['Frequenza originale']
                FREQ_DEG = dictionary[key][sub_key]['Frequenza non sinonimo']
                N_DEGS = dictionary[key][sub_key]['Numero non sinonimi']
                N_TOT = dictionary[key][sub_key]['Numero totali']
                
                Score.append(SCORE)
                Freq_wt.append(FREQ_WT)
                Freq_deg.append(FREQ_DEG)
                N_degs.append(N_DEGS)
                N_tot.append(N_TOT)
                
            except:
                continue
                      
    
    return Score, Freq_wt, Freq_deg, N_degs, N_tot









def dataframe_to_list(dataframe):
    lista = list(dataframe.columns)
    colonne_estratte = []
    if lista[0] == 'Unnamed: 0':
        del lista[0]
    for col in lista:
        temp_col = list(dataframe[col])
        colonne_estratte.append(temp_col)

    return colonne_estratte



def convert_residues3to1(list_of_residues):
    list_converted = []
    for i in list_of_residues:
        y = seq3(i)
        list_converted.append(y)
        
    return list_converted





def create_codon_table():
    bases = "TCAG"
    Codoni = [a + b + c for a in bases for b in bases for c in bases]
    codon_table = dict.fromkeys(Codoni, 0)

    return (codon_table)





    
def count_codons(alignment, Codon_dictionary):
    Lines = alignment.readlines()
    for line in Lines:
        if line[0] == '>':
            continue
        
        else:
            Codons_list = [line[i-3:i] for i in range(0,len(line))]
            for j in Codons_list:
                try:
                    Codon_dictionary[j] += 1
                except:
                    continue

    factor=1.0/sum(Codon_dictionary.values())
    for k in Codon_dictionary:
        Codon_dictionary[k] = Codon_dictionary[k]*factor


    return (Codon_dictionary)





if __name__=='__main__':
    Mutations_file = open(r'/home/blasfemus/Desktop/Nuova_cartella/Lavoro/Synonymous/Working/mavedb.txt')
    Sequence = open(r'/home/blasfemus/Desktop/Nuova_cartella/Lavoro/Synonymous/Working/sequence.txt')
    Degeneri = pd.read_csv(r'/home/blasfemus/Desktop/Nuova_cartella/Lavoro/Synonymous/Working/Tabella_degeneri.csv')
    NS_file = pd.read_csv(r'/home/blasfemus/Desktop/Nuova_cartella/Lavoro/Synonymous/Working/Tabella_non_sinonimi.csv')
    
    
    Pos, Mut, Residues, Scores, PosAA, NS_Pos, NS_Mut, NS_Residues, NS_Scores, NS_PosAA = extract_mutations(Mutations_file)
    RNA_seq, DNA_seq = extract_sequence(Sequence)
    
    Mut_codons, Freq_wt, Freq_mut, Freq_organism_wt, Freq_organism_mut, WT_codons = extract_statistics(RNA_seq, DNA_seq, Mut, Pos)
    
    NS_Mut_codons, NS_Freq_wt, NS_Freq_mut, NS_Freq_organism_wt, NS_Freq_organism_mut, NS_WT_codons = extract_statistics(RNA_seq, DNA_seq, NS_Mut, NS_Pos)
    
    #Mut_codons, WT_codons = extract_statistics(RNA_seq, DNA_seq, Mut, Pos)
    Table_and_Corr(PosAA, Residues, WT_codons, Mut_codons, Freq_wt, Freq_mut, Freq_organism_wt, Freq_organism_mut, Scores, Degeneri, NS_PosAA, NS_Residues, NS_WT_codons, NS_Mut_codons, NS_Scores, NS_file)
    #Alignment = open(r'C:\Users\Andrea\Documents\Lavoro\Synonymous\alignment.fasta')
    #CT = create_codon_table()
    #Dizionario_dei_codoni = count_codons(Alignment, CT)
    #print (Dizionario_dei_codoni)



#factor=1.0/sum(codon_table.values())
#for k in codon_table:
#  codon_table[k] = codon_table[k]*factor





        
