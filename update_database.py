#!/usr/bin/env python
import sys, re, os, getopt
import hashlib
import base64
import sqlite3
import datetime
import happyfile
import time
import progress

count_total_seqs = 0
count_insert_peptide = 0
count_insert_organism = 0
count_insert_taxnode = 0
count_insert_peporg = 0
count_insert_sequence = 0

def xprint(s):
	sys.stderr.write(str(s) + '\n')
	
def xencode(m):
	if sys.version_info < (3,):
		return base64.b64encode(m.digest()).rstrip('=')
	else:
		return base64.encodebytes(m.digest()).decode().replace('\n', '').rstrip('=') 

def seguid(s):
	m = hashlib.sha1()
	m.update(str.encode(s.upper()))
	return xencode(m)

def db_check_create(dbc, type, tbl_name, create_sql):
	dbc.execute("SELECT name FROM sqlite_master WHERE type=? AND name=?", (type, tbl_name))
	if dbc.fetchone():
		return 0
	else:
		dbc.execute(create_sql)
	return 1
	
def create_db_tables(dbc):
	c = 0
	c += db_check_create(dbc, 'table', 'peptide', "CREATE TABLE peptide(pep_id INTEGER PRIMARY KEY, seguid text, seq_id text)")
	if db_check_create(dbc, 'table', 'tax_node', "CREATE TABLE tax_node(tax_id INTEGER PRIMARY KEY, parent_tax_id int, rank int, name text, weight real)") > 0:
		dbc.execute("INSERT INTO tax_node VALUES(0,NULL,0,'root',1)")
		c += 1
	c += db_check_create(dbc, 'table', 'organism', "CREATE TABLE organism(org_id INTEGER PRIMARY KEY, name text, ext_id text, tax_node_id int, tax_str text)")
	c += db_check_create(dbc, 'table', 'pep_org', "CREATE TABLE pep_org(pep_id INTEGER NOT NULL, org_id INTEGER NOT NULL, PRIMARY KEY (pep_id, org_id))")
	c += db_check_create(dbc, 'table', 'sequence', "CREATE TABLE sequence(seq_id text, seq text)")
	db_check_create(dbc, 'index', 'seq_id_index', "CREATE INDEX seq_id_index on peptide(seq_id)")
	db_check_create(dbc, 'index', 'seguid_index', "CREATE INDEX seguid_index on peptide(seguid)")
	if c > 0:
		xprint("[create ] DB tables " + str(c))
		
def store_seq(dbc, seq_id, ext_tax_id, tax_str, seq_str, is_store_seqs):
	global count_insert_peptide
	global count_insert_organism
	global count_insert_peporg
	global count_insert_sequence
	
	if seq_id and ext_tax_id and tax_str and seq_str:
		org_name = tax_str.split(';')[-1]

		# tax_node
		tax_list = store_tax_str(dbc, tax_str)

		tax_list_str = ";".join(str(x) for x in tax_list)
		tax_id = tax_list[-1]
		
		# organism
		dbc.execute("SELECT org_id FROM organism WHERE name = ? AND ext_id = ?", (org_name, ext_tax_id))
		row = dbc.fetchone()
		if row:
			org_id = row[0]
		else:
			dbc.execute("INSERT INTO organism (name, ext_id, tax_node_id, tax_str) VALUES(?, ?, ?, ?)", (org_name, ext_tax_id, tax_id, tax_list_str))
			org_id = dbc.lastrowid
			count_insert_organism += 1
		
		# peptide
		dbc.execute("SELECT pep_id FROM peptide WHERE seguid = ? AND seq_id = ?", (seguid(seq_str), seq_id))
		row = dbc.fetchone()
		if row:
			pep_id = row[0]
		else:
			dbc.execute("INSERT INTO peptide (seguid, seq_id) VALUES(?, ?)", (seguid(seq_str), seq_id))
			pep_id = dbc.lastrowid
			count_insert_peptide += 1
		
		# pep_org
		dbc.execute("SELECT pep_id FROM pep_org WHERE pep_id = ? AND org_id = ?", (pep_id, org_id))
		row = dbc.fetchone()
		if not row:
			dbc.execute("INSERT INTO pep_org (pep_id, org_id) VALUES(?, ?)", (pep_id, org_id))
			count_insert_peporg += 1

		# sequence	
		if is_store_seqs:
			dbc.execute("SELECT seq_id FROM sequence WHERE seq_id = ?", (seq_id,))
			row = dbc.fetchone()
			if not row:
				dbc.execute("INSERT INTO sequence (seq_id, seq) VALUES(?, ?)", (seq_id, seq_str))
				count_insert_sequence += 1

def store_tax_str(dbc, tax_str):
	global count_insert_taxnode
	curr_rank = 0
	parent_id = 0
	tax_node_list = []
	
	for curr_name in tax_str.split(';'):
		curr_rank += 1
		dbc.execute("SELECT tax_id, parent_tax_id FROM tax_node WHERE name = ? AND rank = ?", (curr_name, curr_rank))
		row = dbc.fetchone()
		if row:
			tax_id = row[0]
			parent_id = row[1]
		else:
			dbc.execute("INSERT INTO tax_node (parent_tax_id, rank, name, weight) VALUES(?, ?, ?, ?)", (parent_id, curr_rank, curr_name, 1.0))
			tax_id = dbc.lastrowid
			count_insert_taxnode += 1
		
		tax_node_list.append(tax_id)
		parent_id = tax_id

	return tax_node_list
	
def test_db():
	try:
		dbconn = sqlite3.connect('test.db')
		dbc = dbconn.cursor()
	except:
		xprint("[create_database] test_db: failed")
		return False
		
	xprint("[create_database] test_db: passed")
	return True
		
def test_seguid():
	if seguid("XCHGASCHTHASGCJHGJHJASCTASYCJASHHJASJHSAGDJHADATSF") == "HlN4r5GDkBSYh4TcCDqGzb7ZOZ8":
		xprint("[create_database] test_seguid: passed")
		return True
		
	xprint("[create_database] test_seguid: failed")
	return False
	
def test_all():
	if not (test_db() and test_seguid()):
		sys.exit(2)
		
def main(argv):
	global count_total_seqs
	help = "\n".join([
        "update_database v0.1 (Apr 5, 2017)",
        "Create/update SQLITE3 database from FASTA",
        "",
        "Usage: " + os.path.basename(argv[0]) + " (options) [FASTA file]",
        "   -o file        : output database file (default: LPI_data.db)",
        "   -h, --help     : help", 
        "",
		"Input FASTA headers must have space delimited: >seq_id, taxon_id, taxonomy string",
		"Compression supported: gzip (.gz), bzip2 (.bz2)", ""])

	db_file = 'LPI_data.db'
	fasta_file = ""
	is_store_seqs = False

	try:
		opts, args = getopt.getopt(argv[1:], "o:hs", ["help", "seq", "test"])
	except getopt.GetoptError:
		xprint(help)
		sys.exit(2)

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			xprint(help)
			sys.exit()
		elif opt == '--test':
			test_all()
			sys.exit()
		elif opt == '-o':
			db_file = arg
		elif opt in ("-s", "--seq"):
			is_store_seqs = True 

	if len(args) > 0: 
		fasta_file = args[0]
	else:
		xprint(help)
		sys.exit(2)

	xprint(datetime.datetime.now())
	
	dbconn = sqlite3.connect(db_file)
	dbc = dbconn.cursor()
	xprint("[connect] " + db_file)
	
	create_db_tables(dbc)
	
	fa_handle = happyfile.hopen_or_else_any(fasta_file)
	xprint("reading:  " + fasta_file)

	max_bytes = os.stat(fasta_file).st_size
	curr_bytes = 0
	start_time = time.time()
	progress.draw_progress(curr_bytes, max_bytes, start_time)

	id = ""
	tax_id = ""
	tax_str = ""
	seq = ""
	lndraw = 0
	while 1:
		line = fa_handle.readline()
		if not line:
			break
		curr_bytes += len(line)
		line = line.rstrip()
		lndraw += 1
		
		if line.startswith(">"):
			count_total_seqs += 1
			store_seq(dbc, id, tax_id, tax_str, seq, is_store_seqs)
			m = re.search('^>(\S+)\s+(\S+)\s+(.+)$', line)
			if m:
				id = m.group(1)
				tax_id = m.group(2)
				tax_str = m.group(3)
			seq = ""
		else:
			seq += re.sub('\s', '', line)
		
		if lndraw > 1000:
			progress.draw_progress(curr_bytes, max_bytes, start_time)
			lndraw = 0
			
	store_seq(dbc, id, tax_id, tax_str, seq, is_store_seqs)
	fa_handle.close()
	dbconn.commit()
	dbconn.close()
	
	progress.draw_progress(max_bytes, max_bytes, start_time)
	xprint("")
	
	xprint("total seqs:         " + str(count_total_seqs))
	xprint("[insert ] peptide   " + str(count_insert_peptide))
	xprint("[insert ] organism  " + str(count_insert_organism))
	xprint("[insert ] pep_org   " + str(count_insert_peporg))
	xprint("[insert ] tax_node  " + str(count_insert_taxnode))
	if is_store_seqs:
		xprint("[insert ] sequence  " + str(count_insert_sequence))
	xprint(datetime.datetime.now())
	
if __name__ == "__main__":
    main(sys.argv)
