// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// Created at Sept 20, 2018 
// Updated at Dec 19, 2019
// Updated by Xiaoquan Su
// Version 3.4.2 or above with _OTU_Parser
// For MetaPhlAn

#ifndef sp_parser_h
#define sp_parser_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include "hash.h"
#include "db.h"

#define BUFF 10000000
#define STRLEN 1000
#define TAXA_LEVEL 8 //K, P, C, O, F, G, S, SP

using namespace std;

class _SP_Parser{
      
      public:
             
             _SP_Parser(){};
             
             _SP_Parser(_PMDB db){
                              Database = db;
                              Database.Read_Taxonomy(SP_taxa_table);
                              Load_id();
                              Load_taxon_to_sp();
                              };
                              
             float Load_sp_to_hash(hash_map<string, float, std_string_hash> taxon_count, hash_map<string, float, std_string_hash> & sp_count);
            
      private:
              _PMDB Database;
              hash_map<string, string, std_string_hash> SP_taxa_table;
              
              hash_map <string, vector <int>, std_string_hash> Taxon_to_sp [TAXA_LEVEL - 3];  
              hash_set <string, std_string_hash> Id_hash; //Get_Tree_Id() /species
              vector <string> Id;
              
              int Load_id();
              int Load_taxon_to_sp();
              string Get_taxon(string sp);
             
      };
         
float _SP_Parser::Load_sp_to_hash(hash_map<string, float, std_string_hash> taxon_count, hash_map<string, float, std_string_hash> & sp_count){
                 
                 float count = 0;
                                                 
                 for (hash_map<string, float, std_string_hash>::iterator miter = taxon_count.begin(); miter != taxon_count.end(); miter ++){
                                                               
                     if (Id_hash.count(miter->first) != 0){
                                                 
                                                 if (sp_count.count(miter->first) == 0) sp_count[miter->first] = miter->second;
                                                 else sp_count[miter->first] += miter->second;
                                                 
                                                 }
                     else {
                          if ((miter->first).size() <= 3) continue;
                          if ((miter->first).substr(0, 3) != "s__") continue;
                          
                          string taxon = Get_taxon(miter->first); 
                          bool is_match = false;
                          
                          for (int i = TAXA_LEVEL - 4; i >= 0; i --){
                              if (Taxon_to_sp[i].count(taxon) == 0) continue;
                              is_match = true;
                              vector <int> species = Taxon_to_sp[i][taxon];
                              for (int i = 0; i < species.size(); i ++){
                              
                                  if (sp_count.count(Id[species[i]]) == 0) sp_count[Id[species[i]]] = miter->second / species.size();
                                                 else sp_count[Id[species[i]]] += miter->second / species.size();                              
                                                 }
                          }
                          if (!is_match){                          
                                string nc = "s_unclassified";
                                if (sp_count.count(nc) == 0) sp_count[nc] = miter->second;
                                                 else sp_count[nc] += miter->second;
                                                 }
                                                                              
                          }                 
                     count += miter->second;
                     }
                 
                 return count;
                 }
                                                      

int _SP_Parser::Load_id(){
     
     ifstream infile(Database.Get_Tree_Id().c_str(), ifstream::in);
     if (!infile){
                 cerr << "Error: Cannot open file : " << Database.Get_Tree_Id() << endl;
                 return 0;
                 }
     
     string buffer;
     int count = 0;
     while(getline(infile, buffer)){
                           if (buffer.size() == 0) continue;
                           Id_hash.insert(buffer);
                           Id.push_back(buffer);
                           count ++;
                           }
     
     infile.close();
     infile.clear();
     return count;
     }                    

int _SP_Parser::Load_taxon_to_sp(){
     
     string id_taxon_file = Database.Get_Path() + "/tree/id_taxon.txt";
     
     ifstream infile(id_taxon_file.c_str(), ifstream::in);
     if (!infile){
                 cerr << "Error: Cannot open file : " << id_taxon_file << endl;
                 return 0;
                 }
     
     string buffer;
     int count = 0;
     while(getline(infile, buffer)){
                           if (buffer.size() == 0) continue;
                           stringstream strin(buffer);
                           
                           for (int i = 0; i < TAXA_LEVEL - 3; i++){
                               string taxon;
                               strin >> taxon;
                               Taxon_to_sp[i][taxon].push_back(count);
                               }
                           count ++;
                           }
     
     infile.close();
     infile.clear();
     return count;
     }                    

string _SP_Parser::Get_taxon(string sp){
              
       if (sp.size()<= 3) return "unclassified";
       string taxon = sp.substr(3, sp.size()-3);
       
       string nc = "_unclassified";
       int pos = taxon.find(nc);
       if (pos != string::npos)
          taxon = taxon.substr(0, pos);
       return taxon;
       }

#endif
