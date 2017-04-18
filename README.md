# LPIclassify
#### Taxonomic classification using an extended Lineage Probabilty Index (LPI)

### *Under active development. Check back later for a stable version...*

Usage
-----

| File | Description |
|------|-------------|
| LPIclassify | Main executable |
|  |  |
| create_database.pl | Create initial peptide taxonomy database |
| update_database.py | Update an existing database, only adds new data |


NOTE: the create database is written in Perl because it handles database transactions appropriately without the need for external dependencies, such as APSW. The update database utility checks every field before inserting new data and will run much slower than the create utility.


```
LPIclassify v0.1 (Apr 5, 2017)
Taxonomic classification of query peptides based on BLAST m8 output

Usage: LPIclassify -i (options)
   -a       : minimum score weight (default: 0.1)
   -d file  : database file (default: LPI_data.db)
   -h       : show help
   -i file  : blast m8 file (required, use '-' for STDIN)
   -o file  : output file (default: STDOUT)
   -s       : blast subject IDs are seguids (default: no, sequence IDs)

```

Lineage Probabilty Index (LPI)
------------------------------
For the original description of LPI, see: Podell, S and Gaasterland, T (2007). [DarkHorse: A method for genome-wide prediction of horizontal gene transfer](http://genomebiology.com/2007/8/2/R16). Genome Biology 8(2):R16

Calculation of LPI is modified here to include:
1. Logistic proportion of top BLAST hits
2. Database structure weighting of individual taxa

The resulting highest scoring LPI taxonomy is a representation of the central taxonomic signal, rather than simply the BLAST hit in the database with the highest score. The LPI score on [0,1] provides a level of confidence, or taxonomic consistency among the top BLAST hits.

Installation
------------

Install all dependencies and then run the make utility to compile from source:
```
make
```

Create the database of peptides and corresponding taxonomy using phyloDB or other FASTA file:
```
create_database.pl -o LPI_data.db your_FASTA_file
```
The format of the FASTA headers must be a space separated list of peptide id, taxonomy id, taxonomy string. The taxonomy string is a semi-colon separated list of taxonomic levels, and may include spaces. The last level is understood to be the full organism name.

Example:
```
>gi161784247-NC_010109 4472 Eukaryota;Archaeplastida;Streptophyta;Alismatales;Araceae;Lemnoideae;Lemna;Lemna minor
```

Dependencies
------------

* Python (https://www.python.org/downloads/)
* GNU C++ (https://gcc.gnu.org/)
* Perl (https://www.perl.org/get.html)
* Sqlite3 (https://sqlite.org/download.html)
