// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// Updated at Dec 29, 2016
// Updated by Xiaoquan Su

#include <iostream>

#include "MetaDB_merge.h"
#include "version.h"

using namespace std;

string Infilename1;
string Infilename2;
string Outfilename;

string DB_filename; //for dev

void printhelp(){
    
    cout << "Meta-Storms 2 version " << Version << endl;
    cout << "Compatible with Parallel-META 3" << endl;
    cout << "Usage : MetaDB-merge [Options] Value" << endl;
    
    cout << "\t[Input and Output options]" << endl;
    cout << "\t  -1 The 1st database name [Required]" << endl;
    cout << "\t  -2 The 2nd database name [Required]" << endl;
    cout << "\t  -o Merged output database name, default is \"database_merge.mdb\"" << endl;
    
    cout << endl << "\t[Other options]" << endl;
    cout << "\t  -h Help" << endl;
    
    exit(0);
}

int Parse_Para(int argc, char * argv[]){
    
    Outfilename = "database_merge.mdb";
    
    int i = 1;
    
    if (argc ==1)
        printhelp();
    
    while(i<argc){
        if (argv[i][0] != '-') {
            cerr << "Argument # " << i << " Error : Arguments must start with -" << endl;
            exit(0);
        };
        switch(argv[i][1]){
            case '1': Infilename1 = argv[i+1]; break;
            case '2': Infilename2 = argv[i+1]; break;
            case 'o': Outfilename = argv[i+1]; Outfilename += ".mdb"; break;
            case 'h': printhelp(); break;
            default : cerr << "Error: Unrec argument " << argv[i] << endl; printhelp(); break;
        }
        i+=2;
    }
}

int main(int argc, char * argv[]){
    
    Parse_Para(argc, argv);
    
    cout << "Database merge starts" << endl;
    
    cout << Merge_Database(Infilename1.c_str(), Infilename2.c_str(), Outfilename.c_str()) << " samples merged" << endl;
    
    cout << "Database merge finished" << endl;
    
    return 0;
}
