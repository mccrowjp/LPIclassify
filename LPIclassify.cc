#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <unistd.h>
#include <sqlite3.h>
#include <sys/time.h>
#include "progress.h"
using namespace std;

typedef int taxid_t;
typedef map<taxid_t, vector<taxid_t> > taxonomy_t;
typedef map<string, vector<taxid_t> > peptide_t;
typedef map<taxid_t, string> taxstr_t;
typedef map<taxid_t, double> taxweight_t;

inline bool file_exists(const string& filename) {
    ifstream f(filename.c_str());
    return f.good();
}

inline int stringint0(string val) {
    int retval = 0;
    if(isdigit(val[0])) {
        retval = stoi(val);
    }
    return retval;
}

inline double score_weight(double x) {
    // logistic function, k=15, x0=0.3, x = fraction of maximal bit score
    return (1.0 - (1.0 / (1.0 + exp(-15 * ((1.0 - x) - 0.3)) ) ) );
}

double get_clock_ms() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (tv.tv_sec * 1000.0) + (tv.tv_usec / 1000.0);
}

struct LPI {
    double score;
    taxid_t tax_id;
    LPI() { score = 0.0; tax_id = 0; }
};

class taxonomy_cache {
protected:
    sqlite3_stmt *stmt;
    taxonomy_t t;
    taxstr_t t_str;
    taxweight_t t_weight;
public:
    taxonomy_cache(sqlite3 *db);
    vector<taxid_t> get_lineage(taxid_t taxid);
    string print_lineage(taxid_t taxid);
    double get_weight(taxid_t taxid);
};

class peptide_cache {
protected:
    sqlite3_stmt *stmt;
    peptide_t p;
public:
    peptide_cache(sqlite3 *db, bool is_seguids);
    ~peptide_cache();
    vector<taxid_t> get_peptide_taxids(string pep_id);
};

struct blast_hit {
    string subject_id;
    double bit_score;
    blast_hit(string sid, double bs) { subject_id = sid; bit_score = bs; }
};

struct blast_record {
    string query_id;
    vector < blast_hit > blast_hits;

    void add_blast_hit(string sid, double bs);
    LPI get_best_LPI(peptide_cache *pep_c, taxonomy_cache *tax_c, double min_score_weight);
};

taxonomy_cache::taxonomy_cache(sqlite3 *db) {
    sqlite3_prepare_v2(db, "SELECT name, rank, parent_tax_id, weight FROM tax_node WHERE tax_id = ?", -1, &stmt, NULL);
}

vector<taxid_t> taxonomy_cache::get_lineage(taxid_t taxid) {
    if(taxid > 0) {
        if(t.count(taxid) == 0) {
            vector<taxid_t> term_list;
            int curr_taxid = taxid;
            int parent_taxid = 0;
        
            term_list.push_back(curr_taxid);
            do {
                parent_taxid = 0;
                if(sqlite3_bind_int(stmt, 1, curr_taxid) != SQLITE_OK) {
                    vector<taxid_t> empty;
                    return empty;
                }
                if(sqlite3_step(stmt) == SQLITE_ROW) {
                    parent_taxid = sqlite3_column_int(stmt, 2);
                    if(parent_taxid > 0) {
                        term_list.push_back(parent_taxid);
                    }
                    if(t_str.count(curr_taxid) == 0) {
                        t_str[curr_taxid] = string(reinterpret_cast<const char*>(sqlite3_column_text(stmt, 0)));
                    }
                    if(t_weight.count(curr_taxid) == 0) {
                        t_weight[curr_taxid] = sqlite3_column_double(stmt, 3);
                        if(t_weight[curr_taxid] > 1.0) t_weight[curr_taxid] = 1.0;
                    }
                }
                sqlite3_reset(stmt);
                sqlite3_clear_bindings(stmt);
                curr_taxid = parent_taxid;
            
            } while(curr_taxid > 0);
        
            t[taxid] = term_list;
        }
    
        return t[taxid];
        
    } else {
        vector<taxid_t> empty;    
        return empty;
    }
}

string taxonomy_cache::print_lineage(taxid_t taxid) {
    vector<taxid_t> term_list = get_lineage(taxid);
    stringstream ss;
    
    if(term_list.size() > 0) {
        bool first = true;
        for(vector<taxid_t>::iterator it = term_list.end()-1; it >= term_list.begin(); --it) {
            if(!first) ss << ';';
            ss << t_str.at(*it);
            first = false;
        }
        return ss.str();
    }
    return "";
}

double taxonomy_cache::get_weight(taxid_t taxid) {
    taxweight_t::iterator it = t_weight.find(taxid);
    if(it != t_weight.end()) {
        return it->second;
    }
    return 1.0;
}

peptide_cache::peptide_cache(sqlite3 *db, bool is_seguids) {
    if(is_seguids)
        sqlite3_prepare_v2(db, "SELECT c.tax_node_id FROM peptide a, pep_org b, organism c WHERE a.pep_id = b.pep_id AND b.org_id = c.org_id AND a.seguid = ?", -1, &stmt, NULL);
    else
        sqlite3_prepare_v2(db, "SELECT c.tax_node_id FROM peptide a, pep_org b, organism c WHERE a.pep_id = b.pep_id AND b.org_id = c.org_id AND a.seq_id = ?", -1, &stmt, NULL);
}

vector<taxid_t> peptide_cache::get_peptide_taxids(string pep_id) {
    vector<taxid_t> retval;
    if(pep_id.length() > 0) {
        if(p.count(pep_id) == 0) {
            if(sqlite3_bind_text(stmt, 1, pep_id.c_str(), pep_id.length(), SQLITE_STATIC) != SQLITE_OK) {
                return retval;
            }
            while(sqlite3_step(stmt) == SQLITE_ROW) {
                retval.push_back(sqlite3_column_int(stmt, 0));
            }
            sqlite3_reset(stmt);
            sqlite3_clear_bindings(stmt);

            p[pep_id] = retval;
        }
        return p[pep_id];
    }
    return retval;
}

peptide_cache::~peptide_cache() {
    sqlite3_finalize(stmt);
}

void blast_record::add_blast_hit(string sid, double bs) {
    blast_hit bh(sid, bs);
    blast_hits.push_back(bh);
}

LPI blast_record::get_best_LPI(peptide_cache *pep_c, taxonomy_cache *tax_c, double min_score_weight) {
    LPI retval;
    map< taxid_t, double > taxid_count;
    double bh_max_bitscore = 0.0;

    retval.score = 0.0;
    retval.tax_id = -1;    

    // cerr << "calc LPI: " << query_id << " " << blast_hits.size() << endl;

    for(vector<blast_hit>::iterator it = blast_hits.begin(); it < blast_hits.end(); ++it) {
        if(it->bit_score > bh_max_bitscore) bh_max_bitscore = it->bit_score;
    }
    for(vector<blast_hit>::iterator it = blast_hits.begin(); it < blast_hits.end(); ++it) {
        double sw = score_weight(it->bit_score / bh_max_bitscore);
        // cerr << it->subject_id.c_str() << ", score:" << sw << endl;
        if(sw >= min_score_weight) {
            double time_execute_start = get_clock_ms();
            vector<taxid_t> taxid_list = pep_c->get_peptide_taxids(it->subject_id);
            for(vector<taxid_t>::iterator it2 = taxid_list.begin(); it2 < taxid_list.end(); ++it2) {
                double tw = tax_c->get_weight(*it2);
                taxid_count[*it2] += (sw * tw);
                // cerr << *it2 << " " << taxid_count[*it2] << endl;
            }
            //cerr << "pep (ms): " << it->subject_id.c_str() << " = " << get_clock_ms() - time_execute_start << endl;
        }
    }
    // cerr << "ids: " << taxid_count.size() << endl;
    
    if(taxid_count.size() == 1) {
        retval.score = 1.0;
        retval.tax_id = taxid_count.begin()->first;
        
    } else if(taxid_count.size() > 1) {
        map< int, double > lev_count;
        map< int, double > term_count;
        
        // cerr << blast_hits.size() << "," << taxid_count.size() << endl;
        
        double time_execute_start = get_clock_ms();
        for(map<int,double>::iterator it1 = taxid_count.begin(); it1 != taxid_count.end(); ++it1) {
            vector<int> term_list = tax_c->get_lineage(it1->first);
            int lev = 0;
            for(vector<int>::iterator it2 = term_list.end()-1; it2 >= term_list.begin(); --it2) {
                lev_count[lev] += taxid_count[it1->first];
                term_count[*it2] += taxid_count[it1->first];
                // cerr << it1->first << " " << *it2 << ":" << lev << " " << term_count[*it2] << "," << lev_count[lev] << " = " << taxid_count[it1->first] << endl;
                lev++;
            }
        }
        //cerr << "taxstrs (ms): " << taxid_count.size() << " = " << get_clock_ms() - time_execute_start << endl;
        
        for(map<int,double>::iterator it1 = taxid_count.begin(); it1 != taxid_count.end(); ++it1) {
            vector< int > term_list = tax_c->get_lineage(it1->first);
            double prob_sum = 0;
            double denom_sum = 0;
            int lev = 0;
            for(vector<int>::iterator it2 = term_list.end()-1; it2 >= term_list.begin(); --it2) {
                prob_sum += (1.0 * term_count[*it2] / (lev_count[lev] * (lev+1)));
                denom_sum += (1.0 / (lev+1));
                // cerr << it1->first << " " << *it2 << " " << lev << " : " << term_count[*it2] << " " << lev_count[lev] << " = " << prob_sum << "," << denom_sum << endl;
                lev++;
            }
            double lpi = prob_sum / denom_sum;
            if(lpi > retval.score) {
                retval.score = lpi;
                retval.tax_id = it1->first;
            }
        }
    }
    
    return retval;
}

int main(int argc, char *argv[]) {
    sqlite3 *db;
    ifstream blastifs;
    
    string dbfile = "LPI_data.db";
    string blastfile;
    string outfile;
    double min_score_weight = 0.1;
    bool showhelp = false;
    bool is_seguids = false;
    int max_lwait = 100;
    
    int opt;
    while ((opt = getopt(argc,argv,"a:d:hi:o:s")) != EOF)
        switch(opt) {
        	case 'a': min_score_weight = stoi(optarg); break;
        	case 'd': dbfile = optarg; break;
        	case 'h': showhelp = true; break;
        	case 'i': blastfile = optarg; break;
        	case 'o': outfile = optarg; break;
        	case 's': is_seguids = true; break;
        	default: cerr << "unrecognized argument: " << opt << endl;
        }
    
    if(blastfile.size() == 0) showhelp = true;
    
    if(showhelp) {
    	cerr << "LPIclassify v0.1 (Apr 5, 2017)" << endl
    		 << "Taxonomic classification of query peptides based on BLAST m8 output" << endl
    		 << endl
    		 << "Usage: LPIclassify -i (options)" << endl
    		 << "   -a       : minimum score weight (default: 0.1)"  << endl
    		 << "   -d file  : database file (default: LPI_data.db)"  << endl
    		 << "   -h       : show help"  << endl
    		 << "   -i file  : blast m8 file (required, use '-' for STDIN)"  << endl
    		 << "   -o file  : output file (default: STDOUT)"  << endl
    		 << "   -s       : blast subject IDs are seguids (default: no, sequence IDs)" << endl
    		 << endl;
    	return 0;
    }
    
    double time_start = get_clock_ms();
    
    if(file_exists(dbfile)) {
        cerr << "opening database: " << dbfile << endl;
        
        int rc = sqlite3_open(dbfile.c_str(), &db);
        if(rc) {
            cerr << "Unable to open database: " << sqlite3_errmsg(db) << endl;
            sqlite3_close(db);
            return(1);
        }

        cerr << "reading file: " << blastfile << endl;

        int total_bytes = 374; // get system file size here
        int curr_bytes = 0;
        int lwait_count = 0;
        progressbar pb(total_bytes);
        pb.draw();
        
        blastifs.open(blastfile);
        string line;
        string last_qid;
        blast_record br;
        taxonomy_cache tax_c(db);
        peptide_cache pep_c(db, is_seguids);
        
        if(blastifs.is_open()) {
            while(getline(blastifs, line)) {
                int pos_s = 0;
                int pos_e = 0;
                string qid;
                string sid;
                int bitscore = 0;
                
                lwait_count++;
                curr_bytes += line.size();
                
                for(int i=0; pos_e != string::npos; ++i) {
                    pos_e = line.find('\t', pos_s);

                    if(i==0) qid = line.substr(pos_s, pos_e - pos_s);
                    if(i==1) sid = line.substr(pos_s, pos_e - pos_s);
                    if(i==11) bitscore = stringint0(line.substr(pos_s, pos_e - pos_s));
                    
                    pos_s = pos_e + 1;
                }
                
                if(qid != last_qid) {
                    if(br.blast_hits.size() > 0) {
                        LPI br_lpi = br.get_best_LPI(&pep_c, &tax_c, min_score_weight);
                        cout << br.query_id << '\t' << br_lpi.score << '\t' << tax_c.print_lineage(br_lpi.tax_id) << endl;
                    }
                    br.blast_hits.clear();
                    br.query_id = qid;
                }

                if(bitscore > 0 && !qid.empty()) {
                    br.query_id = qid;
                    br.add_blast_hit(sid, bitscore);
                }
                
                if(lwait_count > max_lwait) {
                    lwait_count = 0;
                    pb.update(curr_bytes);
                }
                
                last_qid = qid;
            }
            if(br.blast_hits.size() > 0) {
                LPI br_lpi = br.get_best_LPI(&pep_c, &tax_c, min_score_weight);
                cout << br.query_id << '\t' << br_lpi.score << '\t' << tax_c.print_lineage(br_lpi.tax_id) << endl;
            }
            
            blastifs.close();
            pb.update(total_bytes);

        } else {
            cerr << "Unable to open file: " << blastfile << endl;
            return 1;
        }

        sqlite3_close(db);

    } else {
        cerr << "Unable to find database file: " << dbfile << endl;
        return 1;
    }
    
    cerr << "time (ms): " << get_clock_ms() - time_start << endl;
   
    return 0;
}
