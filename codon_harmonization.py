
"""

Priority if codon has GC content less in Sce.
Example GCA can't be lowered


"""

import re
import pandas as pd

from urllib.request import Request, urlopen
from bs4 import BeautifulSoup

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils import GC


def download_codon_usage(species_id,species_name,aa_table):
	url='https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species='+str(species_id)+'&aa='+str(aa_table)+'&style=N'
	req=Request(url,headers={'User-Agent': 'Mozilla/5.0'})
	page=urlopen(req).read()
	#
	soup = BeautifulSoup(page, "html.parser")
	#
	open(species_name+'.txt','w').write(soup.find('pre').text)

def get_table(fname):
	df=pd.DataFrame(columns='codon|amino acid|fraction|frequency per 1000|number'.split('|'))
	for line in open(fname,'r').read().split('\n'):
		if line=='':
			continue
		#break
		for line_segment in line[:-1].split(')'):
			#codon,aa,fraction,frequency,number=list(filter(None,line_segment.replace('(',' ').split(' ')))
			row=list(filter(None,line_segment.replace('(',' ').split(' ')))
			df.loc[len(df)]=row
	return df

download_codon_usage(354,'Aztobacter_vinelandii',1)
download_codon_usage(4932,'Saccharomyces_cerevisiae',1)

df_x=get_table('Aztobacter_vinelandii.txt')
df_sce=get_table('Saccharomyces_cerevisiae.txt')



record=SeqIO.read('GAF_domain.fasta','fasta')

sequence_original=str(record.seq)


open('gene_review.txt','w').write('\t'.join('codon|GC content|amino acid|readable x|readable sce'.split('|'))+'\n')
for codon in re.findall('...',sequence_original):
	#
	codon=codon.replace('T','U')
	aa=df_x[df_x['codon']==codon]['amino acid'].item()
	# get codon stats in current organism
	output_x=[[round(GC(df_x['codon'][index])/100,2)]+ [df_x[col][index] for col in 'codon|fraction'.split('|')] for index in df_x[df_x['amino acid']==aa].index.tolist()]
	readable_x=' - '.join([','.join([str(y) for y in x]) for x in sorted(output_x)])
	# get codon stats in Saccharomyces cerevisiae
	output_sce=[[round(GC(df_sce['codon'][index])/100,2)]+ [df_sce[col][index] for col in 'codon|fraction'.split('|')] for index in df_sce[df_sce['amino acid']==aa].index.tolist()]
	readable_sce=' - '.join([','.join([str(y) for y in x]) for x in sorted(output_sce)])
	# write row to file
	printable='\t'.join([codon,str(round(GC(codon)/100,2)),aa,readable_x,readable_sce])
	open('gene_review.txt','a').write(printable+'\n')


"""

# after curation
optimized=pd.read_excel('gene_review.xlsx')


before=Seq(''.join(optimized.codon.tolist()))
after=Seq(''.join(optimized.harmony.tolist()))

GC(before)
GC(after)

# matched translation
str(before.translate())==str(after.translate())
"""



