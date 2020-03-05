// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// Updated at Dec 6, 2018
// Updated by Xiaoquan Su
// Parse the Microbiome Novelty Score by Meta-Storms res

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
string File_out="query.out.mns";

int Max_match = 10;
int Skip = 0;
float Base = 0;
int Mode = 0;

float Get_Score(vector <float> scores, int mode){ //0: liner 1: exp 2: ave
    
    if (scores.size() == 0) return 0;
    float score_sum = 0;
    float weight_sum = 0;
    for (int i = 0; i < scores.size(); i ++){
        float f = 1;        
        switch (mode){
               case 1: f = pow (2, (scores.size() - i)); break;
               case 2: f = 1; break;
               case 0:
               default: f = (float) (scores.size() - i); break;
               }
               
        score_sum += scores[i] * f;
        weight_sum += f;
        
        }
    
    return 1.0 - score_sum / weight_sum;
    }

int Gen_Index(const char * infilename, const char * outfilename){
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
    
    outfile << "#ID\tMNS" << endl;
    
    string buffer;
    unsigned int sample_count = 0;
    while(getline(infile, buffer)){
                          if (buffer[0] == '#') continue;
                          bool no_hit = false;
                          string Exter_Inf, sample_id;
                          stringstream strin(buffer);
                          strin >> Exter_Inf >> sample_id;
                          outfile << sample_id;
                       
                          vector <float> scores;

                            string match_id;
                            float match_sim;
                            float corr_sim;
                            int match_count = 0;
                            int unskip_match_count = 0;
                            float unskip_match_sim = 0;
        
                            while(strin >> match_id){
                                if (match_id == "No-Hit"){
                                    
                                    no_hit = true;
                                    break;
                                    }
                                strin >> match_sim;
                                match_count++;
                                corr_sim = match_sim - Base;
                                corr_sim = (corr_sim > 0)? corr_sim : 0;                                                                
          
                                if (corr_sim == 0) continue;
                                if (match_count <= Skip) continue;
                                
                                scores.push_back(corr_sim);
                                unskip_match_count++;
                                unskip_match_sim += corr_sim;
                                
                                if (unskip_match_count > Max_match) break;
                                }
                
                            if (unskip_match_count == 0) no_hit = true;
                            if (unskip_match_sim == 0) no_hit = true;                                                        
        
                            if (!no_hit)
                                outfile << "\t" << Get_Score(scores, Mode) << endl;
                            else                                
                                outfile << "\t" << 0 << endl;                            
        
                        sample_count ++;
                        }
                          
    infile.close();
    infile.clear();
    
    return sample_count;
    }
    
void PrintHelp(){
     cout << "Meta-Storms 2 version " << Version << endl;
     cout << "Usage : MetaDB-parse-mns [Options] Value" << endl;
     cout << "\t[Input and Output options]" << endl;
     cout << "\t  -i Input file name (the output of MetaDB-search) [Required]"<<endl;
     cout << "\t  -o Output file name, default is \"query.out.mns\"" <<endl;
     
     cout << endl << "\t[Advanced options]" << endl;
     cout << "\t  -b Base of the similarity in the input file, default is 0" << endl;
     cout << "\t  -n Max number of matches in the input file, default is 10" << endl;
     cout << "\t  -s Number of skipped matches in the input file, default is 0" << endl;
     //cout << "\t-M (Upper) Scoring mode: 0: Liner, 1: Exp, 2: Mean, default is 0" << endl;
      
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
                                                      
                  case 'b': Base = atof(argv[i+1]); break;
                  case 'n': Max_match = atoi(argv[i+1]); break;
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
    
    cout << "Microbiome Novelty Score parser" << endl;
    
    Gen_Index(File_ndx1.c_str(), File_out.c_str());
        
    }
