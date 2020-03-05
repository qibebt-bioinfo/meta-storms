// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// Updated at Dec 6, 2018
// Updated by Xiaoquan Su
// Parse the Microbiome Attention Score by Meta-Storms res

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <stdlib.h>
#include <math.h>

#include "version.h"

using namespace std;

string File_ndx1;
string Meta_file;

string File_out="query.out.mas";

int Max_match = 10;
int Skip = 0;
float Base = 0;
int Meta_column = 0;

map <string, string> Meta_data_match;

int Load_Meta(const char * infilename, map <string, string> & meta_data, int c){
    ifstream infile(infilename, ifstream::in);
    if (!infile){
        cout<<"Error: Cannot open the control file : " << infilename << endl;
        return 0;
    }
    
    //unsigned int count = 0;
    string buffer;
    getline(infile, buffer); //title
    
    while(getline(infile, buffer)){
        stringstream strin(buffer);
        string id, meta;
        strin >> id;
        int c_count = 0;
        while(strin >> meta){
         if (c_count == c) break;
         c_count ++;
         }
        meta_data[id] = meta;
        }
   
    infile.close();
    infile.clear();
    return meta_data.size();
    }

int Gen_Index(const char * infilename, const char * outfilename){
    
    map <string, int> match;
    map <string, float> mas;
    vector <string> order;
    
    ifstream infile(infilename, ifstream::in);
    if (!infile){
                 cerr << "Error: Cannot open the input file : " << infilename << endl;
                 return 0;
                 }
          
    string buffer;
    unsigned int sample_count = 0;
    while(getline(infile, buffer)){
                          if (buffer[0] == '#') continue;
                          string Exter_Inf, sample_id;
                          stringstream strin(buffer);
                          strin >> Exter_Inf >> sample_id;
                          order.push_back(sample_id);
                          string sample_id_meta = "QNA";
                          //get query meta
                          if (Meta_data_match.count(sample_id) != 0)
                             sample_id_meta = Meta_data_match[sample_id];
                              
                            string match_id;
                            float match_sim;
                            
                            int match_count = 0;   
                            int unskip_match_count = 0;                            
        
                            while(strin >> match_id){
                                if (match_id == "No-Hit"){
                                                                        
                                    break;
                                    }
                                
                                string match_id_meta = "MNA";
                                if (Meta_data_match.count(match_id) != 0)
                                   match_id_meta = Meta_data_match[match_id];
                                   
                                strin >> match_sim;
                                match_count++;        
                                if (match_count <= Skip) continue;
                                if (match_sim < Base) continue;
                                if (match_id_meta == sample_id_meta) continue;
                                
                                unskip_match_count ++;     
                                if (unskip_match_count > Max_match) break;
                                
                                if (match.count(match_id) == 0){
                                                          
                                                          match[match_id] = 1;
                                                          mas[match_id] = match_sim;
                                                          }
                                else{
                                     match[match_id] ++;
                                     mas[match_id] += match_sim;                                                                          
                                     }
                                }
                                                    
                        sample_count ++;
                        }
                          
    infile.close();
    infile.clear();
    
    
    ofstream outfile(outfilename, ofstream::out);
    if (!outfile){
            cerr << "Error: Cannot open the output file : " << outfilename << endl;
            return 0;
            }
    
    outfile << "ID\tMAS" << endl;
    
    for (int i = 0; i < order.size(); i ++){
        
        outfile << order[i] << "\t";

       if (match.count(order[i]) == 0)
           outfile << 0 << endl;
        else
            outfile << mas[order[i]] << endl; 
        
        }
    
    outfile.close();
    outfile.clear();
    
    return sample_count;
    }
    
void PrintHelp(){
     
     cout << "Meta-Storms 2 version " << Version << endl;
     cout << "Usage : MetaDB-parse-mas [Options] Value" << endl;
     cout << "Basic Options : " << endl;
     cout << "Options :" <<endl;
     
     cout << "\t[Input and Output options]" << endl;
     cout << "\t  -i Input file name (the output of MetaDB-search) [Required]"<<endl;
     cout << "\t  -o Output file name, default is \"query.out.mas\"" <<endl;
     
     cout << endl << "\t[Advanced options]" << endl;
     cout << "\t  -b Base of the similarity in the input file, default is 0" << endl;
     cout << "\t  -n Max number of matches in the input file, default is 10" << endl;
     cout << "\t  -s Number of skipped matches in the input file, default is 0" << endl;
     cout << "\t  -m Input Meta-data file name to eliminate dup-bias [Optional]" << endl;
     cout << "\t  -l Meta-data column, default is 1 (exclude the ID column) [for -m]" << endl;
     
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
                  case 'i': File_ndx1  =argv[i+1];break;
           
                  case 'o': File_out   =argv[i+1];break;
                  case 's': Skip = atoi(argv[i+1]); break;
                  case 'b': Base = atof(argv[i+1]); break;
                  case 'n': Max_match = atoi(argv[i+1]); break;
                  case 'm': Meta_file  =argv[i+1]; break;
                  case 'l': Meta_column = atoi(argv[i+1]) - 1; break;
                  
                  case 'h': PrintHelp();
                  default : cerr << "Error: Unrec argument " << argv[i] << endl; PrintHelp(); break; 
                  }
                  i+=2;
          }
          return ;
}

int main(int argc, char *argv[]){
    
    param(argc, argv);
    
    cout << "Microbiome Attention Score parser" << endl;
    if (Meta_file.size() > 0)
       Load_Meta(Meta_file.c_str(), Meta_data_match, Meta_column);
       
    Gen_Index(File_ndx1.c_str(), File_out.c_str());
    }
