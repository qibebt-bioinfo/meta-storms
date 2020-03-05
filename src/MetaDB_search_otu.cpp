// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// Updated at Dec 6, 2018
// Updated by Xiaoquan Su
// Added unweighted
// Biom format

#include <iostream>
#include <unistd.h>
#include "MetaDB_search.h"

using namespace std;

string DB_filename;
string DB_prefix; //For HDD model only
string Singlefilename;
string Q_Listfilename;
string Listprefix;
string Q_Tablefilename;
string Outfilename;

int Mode = 0; //0:singlefile; 1: listfile; 2: tablefile 3: biomfile

bool Is_index;

int Query_N;
int Query_f;

float Index_t;
float Sim_t; 

int Level;
int Out_format;
int Coren;

int Dist_metric;
int Mem_mode;

void printhelp(){
     
     cout << "Meta-Storms 2 version " << Version << endl;
     cout << "Usage : MetaDB-search-otu [Options] Value" << endl;
     cout << "\t[Database options]" << endl;
     cout << "\t  -d Database file (*.mdb) [Required]" << endl;
     cout << "\t  -H (upper) If enable the HDD mode (low RAM usage), T(rue) or F(alse), default is F" << endl;
     cout << "\t  -P (upper) Path for the HDD mode data files [Optional for '-H T']" << endl;
     
     cout << endl << "\t[Input options]" << endl;
     cout << "\t  -i Single input file name" << endl;
     cout << "\tor" << endl;
     cout << "\t  -l Input filename list" << endl;
     cout << "\t  -p Input List file path prefix for '-l' [Optional for -l]" << endl;
     cout << "\tor" << endl;
     cout << "\t  -T (upper) Input OTU table (*.Count)" << endl;
     
     cout << endl <<  "\t[Output options]" << endl;
     cout << "\t  -o Output file, default is \"query.out\"" << endl;
     
     cout << endl << "\t[Advanced options]" << endl;
     cout << "\t  -n Number of the matched sample(s), default is 10" << endl;
     cout << "\t  -m Minimum similarity of the matched sample(s), range (0.0 ~ 1.0], default is 0" << endl;
     cout << "\t  -e If enable the exhaustive search (low speed), T(rue) or F(alse), default is F" << endl;
     cout << "\t  -w Abundance weighted or unweighted, T(rue) or F(alse), default is T" << endl;
     
     cout << endl << "\t[Other options]" << endl;
     cout << "\t  -t CPU core number, default is auto" << endl;
     cout << "\t  -h Help" << endl;
     exit(0);
     }

int Parse_Para(int argc, char * argv[]){
    
    Is_index = true;
    
    Mode = 0;
    
    Query_N = 10;
    Query_f = 20;
    
    Index_t = 0.1;
    Sim_t = 0;
    
    Level = DEFAULT_INDEX_LEVEL;
    Out_format = 0;
    Coren = 0;
    
    Dist_metric = 0;
    Mem_mode = 0;
    
    Outfilename = "query.out";
    
    int i = 1;
      
      if (argc ==1) 
		printhelp();
      
      while(i<argc){
         if (argv[i][0] != '-') {
                           cerr << "Argument # " << i << " Error : Arguments must start with -" << endl;
                           exit(0);
                           };           
         switch(argv[i][1]){
                            case 'i': Singlefilename = argv[i+1]; Mode = 0; break;
                            case 'l': Q_Listfilename = argv[i+1]; Mode = 1; break;
                            case 'p': Listprefix = argv[i+1]; break;
                            case 'T': Q_Tablefilename = argv[i+1]; Mode = 2; break;
                            case 'd': DB_filename = argv[i+1]; break;                            
                            
                            case 'o': Outfilename = argv[i+1]; break;
                            
                            case 'H': if ((argv[i+1][0] == 't') || (argv[i+1][0] == 'T')) Mem_mode = 1; break;
                            case 'P': DB_prefix = argv[i+1]; break;
                            case 'e': if ((argv[i+1][0] == 't') || (argv[i+1][0] == 'T')) Is_index = false; break;
                            case 'n': Query_N = atoi(argv[i+1]); break;
                            case 'm': Sim_t = atof(argv[i+1]); break;
                            case 'w': if ((argv[i+1][0] == 'f') || (argv[i+1][0] == 'F')) Dist_metric = 1; break;
                            
                            case 't': Coren = atoi(argv[i+1]); break;
                            case 'h': printhelp(); break;
                            default : cerr << "Error: Unrec argument " << argv[i] << endl; printhelp(); break; 
                            }
         i+=2;
         }
  
    if (DB_filename.size() < 0){
                           cerr << "Error: Please assign the database (.mdb) by -d" << endl;
                           exit(0);
                           }
    
    int Max_core_number = sysconf(_SC_NPROCESSORS_CONF);
    
    if ((Coren <= 0) || (Coren > Max_core_number)){
                    //cerr << "Core number must be larger than 0, change to automatic mode" << endl;
                    Coren = Max_core_number;
                    }
    
    if ((Sim_t >= 1.0) || (Sim_t < 0)) 
       Sim_t = 0;
  
  }

void Output_Res(const char * outfilename, vector <string> ids, vector <vector <_MetaDB_Res> > res, int out_format){
    
    //format 0 (default): 1 line mode;
    //format 1: multi-line with title;
    //format 2: multi-line without title;
    
    ofstream outfile(outfilename, ios::out);
    if(!outfile){
                 cerr << "Error: Cannot open output file: " << outfilename << endl;
                 return;
                 }
    
    switch (out_format){
        
        case 1:
            for (int i = 0; i < ids.size(); i ++){
                outfile << "#Query\tMatch\tSimilarity\tQuery_OTU_Count\tMatch_OTU_Count" << endl;
                if (res[i].size() == 0)                    
                    outfile << ids[i] << "\tNo-Hit\tNA\tNA\tNA\tNA" << endl;
                    
                else for (int j = 0; j < res[i].size(); j ++)                        
                        outfile << ids[i] << "\t" << res[i][j].Get_Res_Sample_ID() << "\t" << res[i][j].Get_Res_Sim() << "\t" << res[i][j].Get_Query_OTU_Count() << "\t" << res[i][j].Get_Res_OTU_Count()<< endl;
                        
            }
            break;
        
        case 2:
            for (int i = 0; i < ids.size(); i ++){
                if (res[i].size() == 0)
                    outfile << ids[i] << "\tNo-Hit\tNA\tNA\tNA\tNA" << endl;
                
                else for (int j = 0; j < res[i].size(); j ++)                
                     outfile << ids[i] << "\t" << res[i][j].Get_Res_Sample_ID() << "\t" << res[i][j].Get_Res_Sim() << "\t" << res[i][j].Get_Query_OTU_Count() << "\t" << res[i][j].Get_Res_OTU_Count()<< endl;
                
            }

            break;
            
        case 0:
        default:
            outfile << "#Query\tMatch(es) and similarity(s)" << endl;    
            for (int i = 0; i < ids.size(); i ++){
                outfile << "Query:\t" << ids[i];
                if (res[i].size() == 0) outfile << "\tNo-Hit";
                
                else for (int j = 0; j < res[i].size(); j ++)
                    outfile << "\t" << res[i][j].Get_Res_Sample_ID() << "\t" << res[i][j].Get_Res_Sim();
                outfile << endl;
            }
            break;
        }
    
    outfile.close();
    outfile.clear();
    }


int main(int argc, char * argv[]){
    
    Parse_Para(argc, argv);
    
    cout << "Database search starts" << endl;
    cout << "HDD mode: ";
    if (Mem_mode == 0) cout << "disabled" << endl;
    else cout << "enabled" << endl;
    cout << "Core number: " << Coren << endl;
    
    _MetaDB_search metadb_searcher('G');
    
    if (Mem_mode == 1){
       	if (DB_prefix.size() <= 0) DB_prefix = DB_filename + ".hdd";
		metadb_searcher.Set_Sample_Hdd_Prefix(DB_prefix);
   	    if (!Check_Path(DB_prefix.c_str())){
   	    	cerr << "Error: Cannot open hdd data path: " << DB_prefix << " for HDD mode" << endl;
   	    	cerr << "Please make the hdd data by MetaDB-make-otu, or disable the HDD mode (-H F)" << endl;
   	    	return 0;
		   }
		}
		 
	metadb_searcher.Load_Database(DB_filename.c_str());
    
    cout << "Database loaded" << endl;
    cout << metadb_searcher.Get_Sample_Number() << " samples in the database" << endl;
    
    vector <string> queryids;
    vector <vector <_MetaDB_Res> > res;
    
    switch (Mode){
           
    case 0: //Single
            {
                queryids.push_back("Q_sample");
                if (Is_index)
                    res.push_back(metadb_searcher.Index_Search(Singlefilename.c_str(), Query_N, Query_f, Index_t, Sim_t, Dist_metric, Coren));
                else
                    res.push_back(metadb_searcher.Search(Singlefilename.c_str(), Query_N, Sim_t, Dist_metric, Coren));
            }
            break;
        
    case 1: //List
            {    
                vector <string> queryfiles;
                Load_List(Q_Listfilename.c_str(), queryfiles, queryids, Listprefix);
                if (Is_index)
                    res = metadb_searcher.Index_Search(queryfiles, Query_N, Query_f, Index_t, Sim_t, Dist_metric, Coren);
                else
                    res = metadb_searcher.Search(queryfiles, Query_N, Sim_t, Dist_metric, Coren);
            }
            break;
    
    case 2: //Tableformat
            {
                _Table_Format table(Q_Tablefilename.c_str());
                queryids = table.Get_Sample_Names();
                if (Is_index)
                    res = metadb_searcher.Index_Search(&table, Query_N, Query_f, Index_t, Sim_t, Dist_metric, Coren);
                else
                    res = metadb_searcher.Search(&table, Query_N, Sim_t, Dist_metric, Coren);
            }
            break;
    
     default: cerr << "Error: Input format error" << endl;
              break;
    }
    
    Output_Res(Outfilename.c_str(), queryids, res, Out_format);

    cout << "Database search finished" << endl;
    
    return 0;
    }
