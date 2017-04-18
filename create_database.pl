#!/usr/bin/env perl
use strict;
use DBI;
use Digest::SHA1;
use Time::HiRes qw(gettimeofday);
use progress;
use Getopt::Long;

my $fafile;
my $dbfile = "LPI_data.db";
my $block_size = 10000;
my $is_store_seqs = 0;
my $showhelp = 0;
my $batchmode = 0;

my $count_total = 0;
my $count_pep = 0;
my $count_org = 0;
my $count_tax = 0;
my $count_peporg = 0;
my $count_seq = 0;
my $count_skipped = 0;

# db statement handles, set after tables created
my $dbh;
my $sth_get_taxnode;
my $sth_ins_taxnode;
my $sth_get_organism;
my $sth_ins_organism;
my $sth_get_peptide;
my $sth_ins_peptide;
my $sth_get_peporg;
my $sth_ins_peporg;
my $sth_get_sequence;
my $sth_ins_sequence;

my $max_pep_id = 0;
my $max_tax_id = 0;
my $max_org_id = 0;

my %tax_node_hash;
my %tax_parent_hash;
my %org_hash;

###

sub handle_error {
    my $error = shift;
    print STDERR "SQLite error: $error\n";
    return 1;
}

sub get_seguid {
    my $seq = shift;
    $seq =~ tr/a-z/A-Z/;
    my $sha1 = Digest::SHA1->new;
    $sha1->add($seq);
    return $sha1->b64digest;
}

sub create_db_tables {
    print STDERR "[create ] DB tables\n";
    $dbh->do("CREATE TABLE peptide(pep_id INTEGER PRIMARY KEY, seguid text, seq_id text)");
	$dbh->do("CREATE TABLE tax_node(tax_id INTEGER PRIMARY KEY, parent_tax_id int, rank int, name text, weight real)");
	$dbh->do("INSERT INTO tax_node VALUES(0,NULL,0,'root',1)");
	$dbh->do("CREATE TABLE organism(org_id INTEGER PRIMARY KEY, name text, ext_id text, tax_node_id int, tax_str text)");
	$dbh->do("CREATE TABLE pep_org(pep_id INTEGER NOT NULL, org_id INTEGER NOT NULL, PRIMARY KEY (pep_id, org_id))");
	$dbh->do("CREATE TABLE sequence(seq_id text, seq text)");
}

sub create_db_indexes {
	print STDERR "[create ] DB indexes\n";
    $dbh->do("CREATE INDEX seq_id_index on peptide(seq_id)");
	$dbh->do("CREATE INDEX seguid_index on peptide(seguid)");
}

sub store_tax_str {
    my $tax_str = shift;
    my @taxnodelist = ();
    my $rank = 0;
    my $tax_id = 0;
    my $parent_id = 0;
    
    foreach my $t (split(/;/, $tax_str)) {
        $rank++;
        my $tid = $tax_node_hash{$t}{$rank};
        if($tid) {
            $tax_id = $tid;
            $parent_id = $tax_parent_hash{$tax_id};
        } else {
            $max_tax_id++;
            $tax_id = $max_tax_id;
            $sth_ins_taxnode->execute($tax_id, $parent_id, $rank, $t, 1.0);
            $tax_node_hash{$t}{$rank} = $tax_id;
            $tax_parent_hash{$tax_id} = $parent_id;
            $count_tax++;
        }   
        
        push(@taxnodelist, $tax_id);
        $parent_id = $tax_id;  
    }
    return @taxnodelist;
}

sub store_seq {
    my ($seq_id, $ext_tax_id, $tax_str, $seq_str) = @_;
    if($seq_id && $ext_tax_id && $tax_str && $seq_str) {
        my $seguid = get_seguid($seq_str);
        my @taxstr_list = split(/;/, $tax_str);
        my $org_name = $taxstr_list[scalar(@taxstr_list)-1];
        
        my @taxid_list = store_tax_str($tax_str);
        
        my $tax_list_str = join(";", @taxid_list);
        my $tax_id = @taxid_list[scalar(@taxid_list)-1];
        
        # organism
        my $org_id = $org_hash{$org_name}{$ext_tax_id};
        unless($org_id) {
            $max_org_id++;
            $org_id = $max_org_id;
            $sth_ins_organism->execute($org_id, $org_name, $ext_tax_id, $tax_id, $tax_list_str);
            $org_hash{$org_name}{$ext_tax_id} = $org_id;
            $count_org++;
        }
        
        # peptide
        $max_pep_id++;
        my $pep_id = $max_pep_id;
        $sth_ins_peptide->execute($pep_id, $seguid, $seq_id);
        $count_pep++;

        # pep_org
        $sth_ins_peporg->execute($pep_id, $org_id);
        $count_peporg++;
        
        # sequence
        if($is_store_seqs) {
            $sth_ins_sequence->execute($seq_id, $seq_str);
            $count_seq++;
        }
    }
}

###

GetOptions ("b" => \$batchmode,
            "h" => \$showhelp,
            "o=s" => \$dbfile);

my $help = <<HELP;
create_database v0.1 (Apr 5, 2017)		
Create SQLITE3 database from FASTA

Usage: $0 (options) [FASTA file]
    -b      : batch mode (not interactive)
    -o file : output database file (default: LPI_data.db)
    -h      : show help

Input FASTA headers must have space delimited: >seq_id, taxon_id, taxonomy string
Compression supported: gzip (.gz), bzip2 (.bz2)

HELP

$fafile = shift;

if($showhelp || !$fafile) {
    die $help;
}

if(-e $dbfile) {
    my $yn = "";
    print STDERR "Database exists: $dbfile\n";
    if(!$batchmode) {
        print STDERR "Do you want to overwrite the database (y/N)?\n";
        $yn = <STDIN>;
    } else {
        print STDERR "overwriting file (-b is set)\n";
    }
    if($batchmode || $yn =~ /^y/i) {
        system("rm $dbfile");
    } else {
        die "No changes made. Use update_database.py to perform incremental updates.\n";
    }
}

$dbh = DBI->connect("dbi:SQLite:dbname=".$dbfile, "", "",                         
    { RaiseError => 1, HandleError=>\&handle_error }, ) or die $DBI::errstr;

print STDERR localtime()."\n";
open(IN, $fafile) or die "Unable to open file $fafile\n";

my $max_bytes = -s $fafile;
my $curr_bytes = 0;
my $start_time = gettimeofday();

$dbh->begin_work();
create_db_tables();
$dbh->commit();

if(!$batchmode) {
    drawprogress($curr_bytes, $max_bytes, $start_time);
}

$sth_ins_taxnode = $dbh->prepare("INSERT INTO tax_node (tax_id, parent_tax_id, rank, name, weight) VALUES(?, ?, ?, ?, ?)");
$sth_ins_organism = $dbh->prepare("INSERT INTO organism (org_id, name, ext_id, tax_node_id, tax_str) VALUES(?, ?, ?, ?, ?)");
$sth_ins_peptide = $dbh->prepare("INSERT INTO peptide (pep_id, seguid, seq_id) VALUES(?, ?, ?)");
$sth_ins_peporg = $dbh->prepare("INSERT INTO pep_org (pep_id, org_id) VALUES(?, ?)");
$sth_ins_sequence = $dbh->prepare("INSERT INTO sequence (seq_id, seq) VALUES(?, ?)");

my $lndraw = 0;
my $id;
my $tax_id;
my $tax_str;
my $seq;

print STDERR "reading $fafile\n";

$dbh->begin_work();
while(<IN>) {
    $lndraw++;
    $curr_bytes += length($_);
    chomp;
    
    if(/^>/) {
        $count_total++;
        store_seq($id, $tax_id, $tax_str, $seq);
        if(/^>(\S+)\s+(\S+)\s+(.+)$/) {
            $id = $1;
            $tax_id = $2;
            $tax_str = $3;
        } else {
            $id = "";
            $tax_id = "";
            $tax_str = "";
            $count_skipped++;
        }
        $seq = "";
        
    } else {
        s/[^a-zA-Z]//g;
        $seq .= $_;
    }
    
    if($lndraw >= $block_size) {
        $dbh->commit();
        $dbh->begin_work();
        if(!$batchmode) {
            drawprogress($curr_bytes, $max_bytes, $start_time);
        }
        $lndraw = 0;
    }
}
store_seq($id, $tax_id, $tax_str, $seq);
$dbh->commit();
close(IN);

if(!$batchmode) {
    drawprogress($max_bytes, $max_bytes, $start_time);
    print STDERR "\n";
}

$dbh->begin_work();
create_db_indexes();
$dbh->commit();

print STDERR "total seqs:         $count_total\n";
if($count_skipped > 0) {
    print STDERR "bad headers:        $count_skipped\n";
}
print STDERR "[insert ] peptide   $count_pep\n";
print STDERR "[insert ] organism  $count_org\n";
print STDERR "[insert ] pep_org   $count_peporg\n";
print STDERR "[insert ] tax_node  $count_tax\n";
if($is_store_seqs) {
    print STDERR "[insert ] sequence  $count_seq\n";
}
print STDERR localtime()."\n";
