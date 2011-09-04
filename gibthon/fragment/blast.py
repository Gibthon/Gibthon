#handle blast queries over the web
from Bio.Blast import NCBIWWW
#parsing blast output
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio import Alphabet
"""
This file performs various nucleotide blast searches, over the internet
"""

#programs is a list of valid programs
programs = [#'blastp',
			'blastn',
			'blastx',
			#'tblastn',
			'tblastx',]

#databases: lists of valid databases
protein_db = [	'nr',
				'refseq',
				'swissprot',
				'pat',
				'month',
				'pdb',
				'env_nr',]

nucleotide_db = [	'nr',
					'refseq_mrna',
					'refseq_genomic',
					'est',
					'est_human',
					'est_mouse',
					'est_others',
					'gss',
					'htgs',
					'pat',
					'pdb',
					'month',
					'alu_repeats',
					'dbsts',
					'chromosome',
					'wgs',
					'env_nt',]

#searches maps the program to the required input alphabet & database type
searches = {	#'blastp' : (Alphabet.ProteinAlphabet, protein_db),
				'blastn' : (Alphabet.NucleotideAlphabet, nucleotide_db),
				'blastx' : (Alphabet.NucleotideAlphabet, protein_db),
				#'tblastn': (Alphabet.ProteinAlphabet, nucleotide_db),
				'tblastx': (Alphabet.NucleotideAlphabet, nucleotide_db),}
				

def net_blast(query_record, program='blastn', database = 'nr' ):
	"""
	net_blast(query_record, program, database = 'nr')
	*Perform a BLAST search over the net using the specified program & database
	*before searching, check that the search alphabet is compatible with the type of search,
	*raise a ValueError if not
	
	ARGUMENTS
	query_record: a SeqRecord object containing the query sequence
	program: the program to use, as per:
		http://www.ncbi.nlm.nih.gov/BLAST/blast_program.shtml
	database: the db to query, as per:
		http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=ProgSelectionGuide#db
	
	"""
	#check whether we have a valid query
	if not isinstance(query_record, SeqRecord):
		raise ValueError(u'Invalid Search Item')
	if len(query_record.seq) < 10:
		raise ValueError(u"Query sequence is too short")
	#check that the program is valid
	program = program.lower()
	if program not in searches:
		raise ValueError(u"Invalid Program '%s'" % program)
	
	#check that the alphabet and db are ok
	(required_alpha, required_dbs) = searches[program]
	if not isinstance(query_record.seq.alphabet, required_alpha):
		raise ValueError(u"Query alphabet for '%s' must be '%s'" % (program, alphabets[program]))
	if not (database in protein_db or database in nucleotide_db):
		raise ValueError(u"Invalid database '%s'" % database) 
	if not database in required_dbs:
		raise ValueError(u"Database '%s' cannot be used with program '%s'" % (database, program))
	
	#Value checking done, time to run the search
	results = NCBIWWW.qblast(program, database, query_record.seq, format_type='XML')
	
	#parse the results
	blast_records = NCBIXML.parse(results)
	
	return blast_records

def print_blast_record(br):
	"""
	Print the blast records, for debugging
	"""
	for alignment in br.alignments:
		print '#####Alignment#####'
		print 'sequence: %s' % alignment.title
		for hsp in alignment.hsps:
			print '  *******HSP***'
			print '  length: %i' % alignment.length
			print '  e value: %s' % hsp.expect
			print '  %s ...' % hsp.query[0:75]
			print '  %s ...' % hsp.match[0:75]
			print '  %s ...' % hsp.sbjct[0:75]
