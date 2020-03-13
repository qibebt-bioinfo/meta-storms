// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// Updated at Mar 12, 2020
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
    cout << "Usage : MetaDB-merge [Options] Value" << endl;
    
    cout << "\t[Input and Output options]" << endl;
    cout << "\t  -1 The 1st database name [Required]" << endl;
    cout << "\t  -2 The 2nd database name [Required]" << endl;
    cout << "\t  -o Merged output database name, default is \"database_merge.mdb*\"" << endl;
    
    cout << endl << "\t[Other options]" << endl;
    cout << "\t  -h Help" << endl;
    
    exit(0);
}

int Parse_Para(int argc, char * argv[]){
            
    int i = 1;
    
    if (argc ==1)
        printhelp();
    
    Outfilename = "database_merge";
    
    while(i<argc){
        if (argv[i][0] != '-') {
            cerr << "Argument # " << i << " Error : Arguments must start with -" << endl;
            exit(0);
        };
        switch(argv[i][1]){
            case '1': Infilename1 = argv[i+1]; break;
            case '2': Infilename2 = argv[i+1]; break;
            case 'o': Outfilename = argv[i+1]; break;
            case 'h': printhelp(); break;
            default : cerr << "Error: Unrec argument " << argv[i] << endl; printhelp(); break;
        }
        i+=2;
    }
    if (Infilename1[Infilename1.size() - 1] != Infilename2[Infilename2.size() - 1]) {
                                       cerr << "Error: Two databases must be in the same type." << endl;
                                       exit(0);
                                       }
    
    switch (Infilename1[Infilename1.size() - 1]){
           case 'b' : Outfilename += ".mdb"; break;
           case 's' : Outfilename += ".mdbs"; break;
           case 'f' : Outfilename += ".mdbf"; break;
           default: Outfilename += ".mdb"; break;
           }
   return 0;   
}

int main(int argc, char * argv[]){
    
    Parse_Para(argc, argv);
    
    cout << "Database merge starts" << endl;
    
    cout << Merge_Database(Infilename1.c_str(), Infilename2.c_str(), Outfilename.c_str()) << " samples merged" << endl;
    
    cout << "Database merge finished" << endl;
    
    return 0;
}
