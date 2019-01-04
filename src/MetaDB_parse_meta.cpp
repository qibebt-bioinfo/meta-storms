// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// Updated at Dec 6, 2018
// Updated by Xiaoquan Su
// Parse the Microbiome meta data by Meta-Storms res

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
string File_out="query.out.meta";

int Meta_column = 0;
int Max_match = 10;
int Skip = 0;
float Base = 0;
int Mode = 0;
int Res_num = 1;

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

int Gen_Index(const char * infilename, const char * outfilename, int mode){
    ifstream infile(infilename, ifstream::in);
    if (!infile){
                 cerr << "Error: Cannot open the input file : " << infilename << endl;
                 return 0;
                 }
   
    ofstream outfile(outfilename, ofstream::out);
    if (!outfile){
            cerr << "Error: Cannot open the output file : " << outfilename << endl;
            return 0;
            }
    
    
    outfile << "#ID\tMeta-data\tLikelihood_score" << endl;
    
    string buffer;
    unsigned int sample_count = 0;
    while(getline(infile, buffer)){
                          if (buffer[0] == '#') continue;
                          bool no_hit = false;
                          string Exter_Inf, sample_id;
                          stringstream strin(buffer);
                          strin >> Exter_Inf >> sample_id;
                          outfile << sample_id;
                                                                                                                              
                            string match_id;
                            float match_sim;                            
                            
                            vector <string> match_metas;
                            vector <float> match_sims;
                            
                            int match_count = 0;
                            int unskip_match_count = 0;                                                                         
                            
                            while(strin >> match_id){
                            
                                    if (match_id == "No-Hit"){                                    
                                       no_hit = true;
                                       break;
                                    }   
                                    
                                    strin >> match_sim;
                                    match_count++;
                                    
                                    float corr_sim = match_sim - Base;    
                                    corr_sim = (corr_sim > 0)? corr_sim : 0;
                                    
                                    if (corr_sim == 0) continue;
                                    if (match_count <= Skip) continue;
                                    if (Meta_data_match.count(match_id) == 0) continue;                                    
                                    
                                    match_metas.push_back(Meta_data_match[match_id]);
                                    match_sims.push_back(corr_sim);
                            
                                    unskip_match_count++;
                                    if (unskip_match_count > Max_match) break;
                            }
                                                                         
                            if (unskip_match_count == 0) no_hit = true;                                                                   
        
                            if (!no_hit){
                                         
                              map <string, float> scores;                
                              float unskip_match_sim = 0;   
                              for (int i = 0; i < match_sims.size(); i ++){
                                  
                                  float f = 1;
                                  switch (mode){
                                         case 1: f = pow (2, (match_sims.size() - i)); break;
                                         case 2: f = 1; break;
                                         case 0:
                                         default: f = (float) (match_sims.size() - i); break;
                                         }
                                  
                                  if (scores.count(match_metas[i]) == 0)
                                     scores[match_metas[i]] = 0;
                                  scores[match_metas[i]] += f * match_sims[i];
                                  unskip_match_sim += f * match_sims[i];
                                  }     
                              //sort
                              vector <string> match_meta_sort;
                              vector <float> match_score_sort;
                              
                              for (map <string, float> ::iterator miter = scores.begin(); miter != scores.end(); miter ++){                                       
                                       match_meta_sort.push_back(miter->first);
                                       match_score_sort.push_back(miter->second);
                                       }
                              for (int i = 0; i < match_meta_sort.size()-1; i ++)
                                  for (int j = i + 1; j < match_meta_sort.size(); j ++)
                                      if (match_score_sort[i] < match_score_sort[j]){
                                                         
                                                         string temp_meta = match_meta_sort[i];
                                                         float temp_score = match_score_sort[i];
                                                         
                                                         match_meta_sort[i] = match_meta_sort[j];
                                                         match_score_sort[i] = match_score_sort[j];
                                                         
                                                         match_meta_sort[j] = temp_meta;
                                                         match_score_sort[j] = temp_score;
                                                         }                              
                              //output
                              
                              int output_num = (Res_num < match_meta_sort.size()) ? Res_num : match_meta_sort.size();
                              
                              for (int i = 0; i < output_num; i ++)
                                  outfile << "\t" << match_meta_sort[i] << "\t" << match_score_sort[i] / unskip_match_sim;
                                outfile << endl;
                            }
                            
                            else                                
                                    outfile << "\tNo-Hit" << endl;                            
        
                        sample_count ++;
                        }
                          
    infile.close();
    infile.clear();
    
    return sample_count;
    }
    
void PrintHelp(){
     cout << "Meta-Storms 2 version " << Version << endl;
     cout << "Compatible with Parallel-META 3" << endl;
     cout << "Usage : MetaDB-parse-meta [Options] Value" << endl;
     cout << "\t[Input and Output options]" << endl;
     cout << "\t  -i Input file name (the output of MetaDB-search) [Required]"<<endl;
     cout << "\t  -m Input meta-data file name (meta-data of the database in MetaDB-search) [Required]" << endl;
     cout << "\t  -l Meta-data column, default is 1 (exclude the ID column)" << endl;
     cout << "\t  -o Output file name, default is \"query.out.meta\"" <<endl;
          
     cout << endl << "\t[Advanced options]" << endl;
     cout << "\t  -r Number of predictd meta-data, default is 1" << endl;     
     cout << "\t  -b Base of the similarity in the input file, default is 0" << endl;
     cout << "\t  -n Max number of matches in the input file, default is 10" << endl;
     cout << "\t  -s Number of skipped matches in the input file, default is 0" << endl;
     //cout << "\t  -M (Upper, Score mode) 0: Liner, 1: Exp, 2: Mean, default is 0" << endl;
     
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
                  case 'm': Meta_file  =argv[i+1]; break;
                  case 'l': Meta_column = atoi(argv[i+1]) - 1; break;
                  case 'o': File_out   =argv[i+1];break;
                  case 'r': Res_num = atoi(argv[i+1]); break;
                  case 'n': Max_match = atoi(argv[i+1]); break;
                  case 'b': Base = atof(argv[i+1]); break;
                  case 's': Skip = atoi(argv[i+1]); break;
                  
                  //case 'M': Mode = atoi(argv[i+1]); break;
                  case 'h': PrintHelp();
                  default : cerr << "Error: Unrec argument " << argv[i] << endl; PrintHelp(); break; 
                  }
                  i+=2;
          }
          return ;
}

int main(int argc, char *argv[]){
    
    param(argc, argv);
    
    cout << "Microbiome Metadata parser" << endl;
    
    Load_Meta(Meta_file.c_str(), Meta_data_match, Meta_column);
    Gen_Index(File_ndx1.c_str(), File_out.c_str(), Mode);
    }
