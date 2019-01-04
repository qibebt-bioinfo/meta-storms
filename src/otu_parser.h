// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// Created at June 9, 2018 
// Updated at Dec 19, 2018
// Updated by Xiaoquan Su
// Version 3.4.2 or above with _OTU_Parser

#ifndef otu_parser_h
#define otu_parser_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include "hash.h"
#include "db.h"

#define BUFF 10000000
#define STRLEN 1000

using namespace std;

class _OTU_Parser{
      
      public:
             
             _OTU_Parser(){};
             
             _OTU_Parser(_PMDB db){
                              Database = db;
                              Database.Read_Taxonomy(OTU_taxa_table);
                              Database.Load_Copy_Number(Cp_number);
                              };
             
             int Output_hash_to_table(const char * tablefilename, hash_map<string, int, std_string_hash> otu_count, bool is_cp){
                 
                 hash_map <string, float, std_string_hash> otu_abd;
                 
                 //calc_abd
                 float seq_count = 0;
                 float otu_count_sum = 0;
                 for (hash_map<string, int, std_string_hash> :: iterator miter = otu_count.begin(); miter != otu_count.end(); miter ++){                     
                                  float cp = 1;
                                  if ((is_cp) && (Cp_number.count(miter->first) != 0))
                                            cp = Cp_number[miter->first];
                                  otu_count_sum += ((float) miter->second / cp );
                                  seq_count += miter->second;
                                  }
                 
                  for (hash_map<string, int, std_string_hash> :: iterator miter = otu_count.begin(); miter != otu_count.end(); miter ++){
                     
                     float cp = 1;
                     if ((is_cp) && (Cp_number.count(miter->first) != 0))
                                 cp = Cp_number[miter->first];
                     
                     otu_abd[miter->first] = ((float) miter->second / cp ) / otu_count_sum;
                     }
                     
                 ofstream outfile(tablefilename, ios::out);
                 if (!outfile){
                               cerr << "Error: Cannot open output file: " << tablefilename << endl;
                               return 0;
                               }
                 
                 outfile << "#Database_OTU\tCount\tAbundance\tTaxonomy" << endl;
                 
                 for (hash_map<string, int, std_string_hash> :: iterator miter = otu_count.begin(); miter != otu_count.end(); miter ++){
                     
                     string taxa = "Unclassified; otu_" + miter->first;
                     if (OTU_taxa_table.count(miter->first) != 0) taxa = OTU_taxa_table[miter->first];
                     
                     outfile << miter->first << "\t" << miter->second << "\t" << otu_abd[miter->first] << "\t" << taxa << endl;
                     }
                                                   
                 outfile.close();
                 outfile.clear(); 
                 
                 return  seq_count;            
                 }
                                           
             int Update_class_taxa(const char * infilename, const char * outfilename){
                 
                 hash_map<string, int, std_string_hash> otu_count;
                 int count = Load_file_to_hash(infilename, otu_count);
                 Output_hash_to_table(outfilename, otu_count, true);
                 
                 return count;   
                 }   
            
            int Load_file_to_hash(const char * classfilename, hash_map<string, int, std_string_hash> & otu_count){
                 
                    //fgets
                    FILE * fptr = fopen(classfilename, "r");
                    if (fptr == NULL){
                       cerr << "Error: Cannot open input file: " << classfilename << endl;
                       return -1;
                       }
                 
                  char * str_buffer = new char [BUFF];
                  vector <string> file_buffer;
                  
                   while(fgets(str_buffer, BUFF-1, fptr) != NULL){
                                        
                                        string a_line = str_buffer;
                                        file_buffer.push_back(a_line);
                                        
                                        }
                   delete [] str_buffer;
                   fclose(fptr);
                 
                 int mode = 1; //default is simple mode
                 
                 string title = file_buffer[0];
                 if (title[1] == 'S') mode = 0;
                 
                 unsigned int count = 0;                              
                 
                 if (mode == 0){ //detail mode, old format                          
                   for(int i = 1; i < file_buffer.size(); i ++){
                   					   char str_seq_id [STRLEN];
                   					   char str_id [STRLEN];
                                       sscanf(file_buffer[i].c_str(), "%s%s", str_seq_id, str_id);
                                       string id = str_id;
                                       
                                       if (otu_count.count(id) == 0)
                                             otu_count[id] = 1;
                                       else otu_count[id] ++;
                                       
                                       count ++;
                                       } 
                     }
                 
                 else{//simple mode, new format 
                       for(int i = 1; i < file_buffer.size(); i ++){
                                       char str_id [STRLEN];
                                       int id_count = 0;                                       
                                       sscanf(file_buffer[i].c_str(), "%s%d", str_id, &id_count);									   
                                       string id = str_id;
                                       otu_count[id] = id_count;
                                       count += id_count;
                                       }                                             
                      } 
                 return count;
                 }
                
      private:

              _PMDB Database;
              hash_map<string, string, std_string_hash> OTU_taxa_table;
              hash_map <string, float, std_string_hash> Cp_number;
      };

#endif
