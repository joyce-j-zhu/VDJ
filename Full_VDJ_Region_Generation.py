import pandas as pd
import numpy as np
import argparse
from difflib import SequenceMatcher
from Bio import Align

#How to use in command line:
#python3 vdj_region_gen_v2.py (input_file) -v (v_file) -j (j_file) -o (output name)

### MAIN FUNCTION ###
def vdj_region_gen(input_file, v_file, j_file, output):
    '''
    Input: vj file, v data file, j data file, output csv file name
    Output: vj file with full V-CDR3-J AA sequences generated
    '''
    #Read in files
    vj_df = pd.read_csv(input_file)
    v_df = pd.read_csv(v_file)
    j_df = pd.read_csv(j_file)
    #Add corresponding V/J AA sequences to the main df
    vj_with_seq = id_match(vj_df, v_df, j_df)
    #Remove rows where AA sequence for VID/JID could not be found
    vj_with_seq_updated = vj_with_seq[(vj_with_seq['V_AA_Seq'] != 'NA') & (vj_with_seq['J_AA_Seq'] != 'NA')]
    #Generate column containing full V-CDR3-J AA sequences
    final_df = full_region_gen(vj_with_seq_updated)
    #Save output to csv
    final_df.to_csv(output, header=True, index=False)

#Function matches VID and JID to extract corresponding AA sequences
def id_match(vj, v, j):
    '''
    Inputs:
        vj = df containing both V and J IDs
        v = df containing V IDs and corresponding AA sequences
        j = df containing J IDs and corresponding AA sequences
    Output:
        df updated with corresponding V and J AA sequences
    '''
    #Convert value to a string to ensure correct type
    vj['vAllele'] = vj['vAllele'].astype(str)
    vj['jAllele'] = vj['jAllele'].astype(str)
    #Create empty lists
    v_seq_list = []
    j_seq_list = []
    #Match V and J of input df to V and J of the V/J files ("v" and "j")
    for i, row in vj.iterrows():
        v_id = row['vAllele']
        j_id = row['jAllele']
        try:
            filtered_v = v[v['V.Name'] == v_id].reset_index()
            v_seq = filtered_v.loc[0]['V.AA.String']
        except:
            v_seq = 'NA'
        try:
            filtered_j = j[j['J.Name'] == j_id].reset_index()
            j_seq = filtered_j.loc[0]['J.AA.String']
        except:
            j_seq = 'NA'
        v_seq_list.append(v_seq)
        j_seq_list.append(j_seq)
    #Add the AA sequences corresponding to Vs and Js to the dataframe
    output_df = vj.copy()
    output_df['V_AA_Seq'] = v_seq_list
    output_df['J_AA_Seq'] = j_seq_list
    return output_df

###Next two functions are used in v_trimming()

#Find all common substrings between two strings
def all_common_substrings(str1, str2):
    common_substrings = set()
    for i in range(len(str1)):
        for j in range(len(str2)):
            k = 0
            while i + k < len(str1) and j + k < len(str2) and str1[i + k] == str2[j + k]:
                k += 1
                common_substrings.add(str1[i:i + k])
    return list(common_substrings)

#Filter substrings that start with 'C' and are at least 2 characters long
def filter_substrings(substrings):
    return [substring for substring in substrings if substring.startswith('C') and len(substring) >= 2]

#Finds longest common sequence between V & CDR3 and keeps only the preceding V sequence
#Requires longest common sequence to start with a 'C'
def v_trimming(v_seq, cdr3_seq):
    '''
    Inputs:
        v_seq = V sequence
        cdr3_seq = CDR3 sequence
    Output:
        v_seq_trimmed = identify lcs between V and CDR3 sequences, locate in V sequence, and trim corresponding end off of V sequence
            1st method: lcs starting with C and with match size of >= 2
            2nd method: pairwise alignment starting with C and with match size of >= 5
    '''
    # Method 1: Longest Common Subsequence (LCS)
    possible_substrings = filter_substrings(all_common_substrings(v_seq[-20:], cdr3_seq))
    if possible_substrings:
        match_seq = max(possible_substrings, key=len)
        v_seq_as_parts = v_seq.split(match_seq)[:-1]
        v_seq_trimmed = match_seq.join(v_seq_as_parts)
    else:
    # Method 2: Pairwise Alignment
        v_segment = v_seq[-20:]
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'

        #Set scoring parameters
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -0.5

        #Perform alignment
        alignments = aligner.align(v_segment, cdr3_seq)

        #Find alignments starting with 'C' and then the one with the greatest match score
        threshold = 0
        filtered_alignments = [
            alignment for alignment in alignments
            if (v_segment[alignment.aligned[0][0][0]].startswith('C'))  #Alignment starts with 'C'
            and (alignment.score >= threshold) #Alignment score meets threshold set above
        ]
        try:
            best_alignment = max(filtered_alignments, key=lambda x: x.score)
            #Extract trimmed V sequence
            new_match_seq = v_segment[best_alignment.aligned[0][0][0]:]
            v_seq_as_parts = v_seq.split(new_match_seq)[:-1]
            v_seq_trimmed = new_match_seq.join(v_seq_as_parts)
        except:
            v_seq_trimmed = 'NA'
    return v_seq_trimmed

#Find motif in J sequence and trim accordingly
def j_trimming(j_seq):
    '''
    Input: j_seq = J sequence
    Output: j_seq_trimmed which is the sequence following F/W motif ((F/W)G(N)G)
    '''
    j_dict = {'F':2, 'W':2, 'G':1}
    value_list = []
    for n in j_seq:
        try:
            value = j_dict[n]
            value_list.append(value)
        except:
            value = 0
            value_list.append(value)
    value_string = ''.join(str(num) for num in value_list)
    match = SequenceMatcher(None, value_string, '2101').find_longest_match()
    if match.size == 4:
        string_start = match.a + 1
        j_seq_trimmed = j_seq[string_start:]
    else:
        match = SequenceMatcher(None, value_string, '2111').find_longest_match()
        if match.size == 4:
            string_start = match.a + 1
            j_seq_trimmed = j_seq[string_start:]
        else:
            j_seq_trimmed = 'NA'
    return j_seq_trimmed

#Applies v_trimming and j_trimming functions to create full region sequence
def full_region_gen(vj):
    '''
    Input: vj file containing V and J AA sequences (output from id_match function)
    Output: df updated with V-CDR3-J full AA sequences
    '''
    #Create empty lists
    v_trim_list = []
    j_trim_list = []
    #Iterate through rows, identify V/CDR3/J sequences, and apply V/J trimming rules
    for i, row in vj.iterrows():
        v = row['V_AA_Seq']
        c = row['CDR3']
        j = row['J_AA_Seq']
        try:
            v_trimmed = v_trimming(v,c)
        except:
            v_trimmed = 'NA'
        try:
            j_trimmed = j_trimming(j)
        except:
            j_trimmed = 'NA'
        v_trim_list.append(v_trimmed)
        j_trim_list.append(j_trimmed)
    #Save trimmed V/J outputs
    vj['V_AA_Trimmed'] = v_trim_list
    vj['J_AA_Trimmed'] = j_trim_list
    #Remove NAs
    vj = vj[(vj.V_AA_Trimmed != 'NA') & (vj.J_AA_Trimmed != 'NA')]
    #Reconstruct full V-CDR3-J sequences
    vj['V_CDR3_J_Sequence'] = vj['V_AA_Trimmed'] + vj['CDR3'] + vj['J_AA_Trimmed']
    return vj

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generates full V-CDR3-J AA sequences')
    parser.add_argument('input_file', type=str, help='CSV file with VID, JID, CDR3 sequence, etc.')
    parser.add_argument('-v', '--v_file', default='V_data.csv', type=str, help='V_data file containing corresponding V AA sequences')
    parser.add_argument('-j', '--j_file', default='J_data.csv', type=str, help='J_data file containing corresponding J AA sequences')
    parser.add_argument('-o', '--output', default='vdj_output.csv', help='Assign name to output CSV file')
    args = parser.parse_args()
    vdj_region_gen(input_file=args.input_file, v_file=args.v_file, j_file=args.j_file, output=args.output)
