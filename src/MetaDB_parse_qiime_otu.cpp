// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// Updated at Dec 7, 2018
// Updated by Xiaoquan Su
// Convert QIIME 1.9.1 otu to PM3 format

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <stdlib.h>
#include "otu_parser.h"

#include "version.h"

using namespace std;

string Infilename;
string Outfilename = "classification.txt";

void PrintHelp(){
     cout << "Meta-Storms 2 version " << Version << endl;
     cout << "Compatible with QIIME 1.9.1" << endl;
     cout << "Usage : MetaDB-parse-qiime-otu [Options] Value" << endl;
     cout << "\t[Input and Output options]" << endl;
     cout << "\t  -i Input file name (the OTU map by QIIME, eg. pick_otu.py) [Required]"<<endl;
     cout << "\t  -o Output file name, default is \"classification.txt\"" <<endl;      
     cout << endl << "\t[Other options]" << endl;
     cout << "\t  -h Help" << endl;
     exit(0);
}

void param(int argc, char *argv[]){
     int i = 1;
     if (argc == 1) PrintHelp();
     while (i<argc){
           if (argv[i][0] != '-'){
              cout << "Argument # " << i;
              cout <<" Error : Arguments must start with -\n"<<endl;
              exit(0);
              }
           switch (argv[i][1]){
                  case 'i': Infilename  =argv[i+1];break;
                  case 'o': Outfilename   =argv[i+1];break;
                  
                  case 'h': PrintHelp();
                  default : cerr << "Error: Unrec argument " << argv[i] << endl; PrintHelp(); break; 
                  }
                  i+=2;
          }
          return ;
}

int Load_qiime_map_to_hash(const char * infilename, hash_map<string, int, std_string_hash> & otu_count){
                 
                 fstream infile(infilename, ios::in);
                 if (!infile){
                              cerr << "Error: Cannot open the input file : " << infilename << endl;
                              return 0;
                              }
                 
                 string buffer;
                 unsigned int line_count = 0;
                 while(getline(infile, buffer)){
                                       if (buffer[0] == '#') continue;
                                       stringstream strin(buffer);
                                       string a_otu;
                                       string seq_id;
                                       strin >> a_otu;
                                       
                                       int a_otu_count = 0;
                                       while(strin >> seq_id)
                                                   a_otu_count ++; 
                                       
                                       if (otu_count.count(a_otu) == 0)
                                                                otu_count[a_otu] = a_otu_count;
                                       else
                                           otu_count[a_otu] += a_otu_count;
                                       
                                       line_count ++;
                                       }
                 
                 infile.close();
                 infile.clear();                    
                 
                 return line_count;
                 }

int main(int argc, char *argv[]){
    
    _PMDB db('g');
    _OTU_Parser otu_parser(db);
    
    hash_map<string, int, std_string_hash> otu_count;
    
    param(argc, argv);
    
    cout << "Parse QIIME OTU" << endl;
    
    Load_qiime_map_to_hash(Infilename.c_str(), otu_count);
    
    cout << otu_parser.Output_hash_to_table(Outfilename.c_str(), otu_count, true) << " sequences are parsed" << endl;
    
    return 0;    
    }
