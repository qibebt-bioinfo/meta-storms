// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// Updated at Oct 25, 2018
// Updated by Xiaoquan Su
// BIOM format

#include <iostream>
#include "MetaDB_make.h"

using namespace std;

string Outfilename;
string Listprefix;
string DB_Listfilename;
string DB_Tablefilename;

int Mode = 0; //0: listfile; 1: tablefile; 2: biomfile 3: only make hdd files 4: dev mode
string DB_filename; //for dev

bool Is_cp_correct;
int Level;

int Mem_mode;
string DB_prefix; //For HDD model only

void printhelp(){
     
    cout << "Meta-Storms 2 version " << Version << endl;
    cout << "Compatible with Parallel-META 3" << endl;
    cout << "Usage : MetaDB-make [Options] Value" << endl;
    
    cout << "\t[Input options]" << endl;
    cout << "\t  -i or -l Input filename list" << endl;
    cout << "\t  -p List file path prefix for '-i' or '-l' [Optional for -i and -l]" << endl;
    cout << "\tor" << endl;
    cout << "\t  -T (upper) Input OTU table (*.Count)" << endl;
    cout << "\tor" << endl;
    cout << "\t  -d (*.mdb) Make the HDD mode data files for a database" << endl;
    
    cout << endl <<  "\t[Output options]" << endl;
    cout << "\t  -o Output database name, default is \"database.mdb\"" << endl;
    cout << "\t  -H (upper) If enable the HDD mode (low RAM usage), T(rue) or F(alse), default is F" << endl;
   
    cout << endl << "\t[Other options]" << endl;
    cout << "\t  -h Help" << endl;
    
    exit(0);      
    }

int Parse_Para(int argc, char * argv[]){
    
    Is_cp_correct = true;
    Level = DEFAULT_INDEX_LEVEL;
    Outfilename = "database.mdb";
    
    Mem_mode = 0;
    
    int i = 1;
      
      if (argc ==1) 
		printhelp();
      
      while(i<argc){
         if (argv[i][0] != '-') {
                           cerr << "Argument # " << i << " Error : Arguments must start with -" << endl;
                           exit(0);
                           };           
         switch(argv[i][1]){
                            case 'l':
                            case 'i': DB_Listfilename = argv[i+1]; Mode = 0; break;   
							case 'p': Listprefix = argv[i+1]; break;   
                            case 'T': DB_Tablefilename = argv[i+1]; Mode = 1; break;    
                                                 
                            case 'd': DB_filename = argv[i+1]; Mode = 3; Mem_mode = 1; break;                                               
                            case 'o': Outfilename = argv[i+1]; Outfilename += ".mdb"; break;
                            //case 'r': if ((argv[i+1][0] == 'f') || (argv[i+1][0] == 'F')) Is_cp_correct = false; break;
                            case 'H': if ((argv[i+1][0] == 't') || (argv[i+1][0] == 'T')) Mem_mode = 1; break;
                            //case 'L': Level = atoi(argv[i+1]); break; //For dev only
                            case 'h': printhelp(); break;
                            default : cerr << "Error: Unrec argument " << argv[i] << endl; printhelp(); break;
                            }
         i+=2;
         }
  }
  
int main(int argc, char * argv[]){
    
    Parse_Para(argc, argv);
    
    cout << "Database construction starts" << endl;
    cout << "HDD mode: ";
    if (Mem_mode == 0) cout << "disabled" << endl;
    else cout << "enabled" << endl;
    
    _MetaDB_make metadb_maker;
    
    if (Mem_mode == 1){
       if (Mode == 3) DB_prefix = DB_filename + ".hdd";
       else DB_prefix = Outfilename + ".hdd";
  	   Check_Path(DB_prefix.c_str(), 1);
	   metadb_maker.Set_Sample_Hdd_Prefix(DB_prefix);
   	   }
   	   
    int count = 0;
    
    switch (Mode){
           case 0:  //List
                    count = metadb_maker.Make_Database(Outfilename.c_str(), DB_Listfilename, Listprefix, Is_cp_correct, Level); 
                    break;
           case 1: //Table
                   {
                   _Table_Format table(DB_Tablefilename.c_str());
                   count = metadb_maker.Make_Database(Outfilename.c_str(), &table, Is_cp_correct, Level);
                   }
                   break;
           case 3: //Hdd only
                   count = metadb_maker.Make_Database_Hdd(DB_filename.c_str());
                   break;
           case 4: //for dev
                   count = metadb_maker.Make_Database(Outfilename.c_str(), DB_filename.c_str(), Level);
                   break;
                
           default: cerr << "Error: Input format error" << endl;
                    break;
           }
	        
    cout << count << " samples loaded" << endl;
    
    cout << "Database construction finished" << endl;
    
    return 0;
    }
