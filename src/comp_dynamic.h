// Updated at July 29, 2019
// Updated by Xiaoquan Su
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// version 3.1 or above with _Table_Format
// version 3.4.1 or above with OTU_Parser
// Added unweighted based on 3.4.1
// For MetaPhlAn

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>

#include "hash.h"
#include "utility.h"
#include "version.h"
#include "table_format.h"
#include "db.h"
#include "otu_parser.h"
#include "sp_parser.h"
#include "comp.h"
#include "dist.h"

#ifndef COMP_DYNAMIC_H
#define COMP_DYNAMIC_H

using namespace std;

class _Comp_Tree_Dynamic : public _Comp_Tree {
      
      public:
             _Comp_Tree_Dynamic() : _Comp_Tree() {};
    
            _Comp_Tree_Dynamic(char db) : _Comp_Tree(db){
                        
                        //Database.Set_DB(db);                                                  
                        Sp_parser = _SP_Parser(Database);                                             
                        }

             int Load_abd(const char * infilename, float * Abd, bool is_cp_correct);
             int Load_abd(const char * infilename, float * Abd);
             int Load_abd(_Table_Format * table, float * Abd, int sample, bool is_cp_correct); //Load by table_format
             int Load_abd(_Table_Format * table, float * Abd, int sample); //Load by table_format

             
      private:
               _SP_Parser Sp_parser;
    
              };

int _Comp_Tree_Dynamic::Load_abd(const char * infilename, float * Abd, bool is_cp_correct){
    
    memset(Abd, 0, LeafN * sizeof(float));
    
    hash_map<string, float, std_string_hash> taxon_count;
    hash_map<string, float, std_string_hash> otu_count;
    hash_map<string, float, std_string_hash> otu_abd;
     
    Otu_parser.Load_file_to_hash(infilename, taxon_count);
    Sp_parser.Load_sp_to_hash(taxon_count, otu_count);
    
    float total = 0;
        
    //cp_number_correct
    
    for (hash_map<string, float, std_string_hash>::iterator miter = otu_count.begin(); miter != otu_count.end(); miter ++){
        
        float cp_no = 1.0;        
        if (is_cp_correct)
           cp_no = Otu_parser.Get_cp_by_OTU(miter->first);
           
        otu_abd[miter->first] = (float) miter->second / cp_no;
        total += otu_abd[miter->first];
        }
        
    //norm
    total /= 100.0;
    int mapped_otu_count = 0;
    for (hash_map<string, float, std_string_hash>::iterator miter = otu_abd.begin(); miter != otu_abd.end(); miter ++){
        miter->second /= total;
        //debug
        //cout << miter->first << endl;
        
        if (Id_hash.count(miter->first) != 0){
            Abd[Id_hash[miter->first]] = miter->second;
            mapped_otu_count ++;
            }
        }
        
    return mapped_otu_count;
    }

int _Comp_Tree_Dynamic::Load_abd(const char * infilename, float * Abd){
    
    return this->Load_abd(infilename, Abd, true);
    }

int _Comp_Tree_Dynamic::Load_abd(_Table_Format * table, float *Abd, int sample, bool is_cp_correct){
    
    memset(Abd, 0, LeafN * sizeof(float));

    vector <string> otus = table->Get_Feature_Names();
    vector <float> abds = table->Get_Abd(sample);
    
    hash_map<string, float, std_string_hash> taxon_count;
    hash_map<string, float, std_string_hash> otu_count;
    hash_map<string, float, std_string_hash> otu_abd;
    
    //load into hash
    for (int i = 0; i < otus.size(); i ++)
        if (abds[i] > 0){
        
           string a_otu = Check_SP(otus[i]);
        
           if (taxon_count.count(a_otu) == 0) taxon_count[a_otu] = abds[i];
           else taxon_count[a_otu] += abds[i];
           }
        
    Sp_parser.Load_sp_to_hash(taxon_count, otu_count); 
    
    float total = 0;
    
    //cp_number_correct
    
    for (hash_map<string, float, std_string_hash>::iterator miter = otu_count.begin(); miter != otu_count.end(); miter ++){
        
        float cp_no = 1.0;        
        if (is_cp_correct)
           cp_no = Otu_parser.Get_cp_by_OTU(miter->first);
           
        otu_abd[miter->first] = (float) miter->second / cp_no;
        total += otu_abd[miter->first];
        }
    //norm
    total /= 100.0;
    int mapped_otu_count = 0;
    for (hash_map<string, float, std_string_hash>::iterator miter = otu_abd.begin(); miter != otu_abd.end(); miter ++){
                          miter->second /= total;
                          if (Id_hash.count(miter->first) != 0){
                             Abd[Id_hash[miter->first]] = miter->second;
                             mapped_otu_count ++;
                             }
                          }
    return mapped_otu_count;
    }

int _Comp_Tree_Dynamic::Load_abd(_Table_Format * table, float *Abd, int sample){
                                       return this->Load_abd(table, Abd, sample, true);
                                       }

#endif
