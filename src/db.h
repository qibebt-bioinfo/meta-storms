// Updated at Dec 20, 2016
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// Code by: Honglei Wang, Gongchao Jing, Xiaoquan Su
// Define the version of Parallel-META toolkit

#include <iostream>
#include "utility.h"
#include "hash.h"

using namespace std;

#ifndef _DB_H
#define _DB_H

class _PMDB{
      
      public:
             _PMDB(){                                         
                     Set_DB('g');                     
                     }
                     
             _PMDB(char db){                                                
                        Set_DB(db);
                        }                        
             
             void Set_DB(char db){
                  
                  DB_path = Check_Env();
                  
                  switch (db){                               
                               case 'e':
                               case 'E': DB_id = 'E';
                                         DB_domain = 1;
                                         DB_name = "sliva_18s"; 
                                         DB_description = "Eukaryote (18S rRNA)"; 
                                         Is_taxa = true; 
                                         Is_tree = true; 
                                         Is_cp = false;
                                         Is_func = false; 
                                         break;
                                         
                               case 'o':
                               case 'O': DB_id = 'O';
                                         DB_domain = 1;
                                         DB_name = "oral_core"; 
                                         DB_description = "Oral Core (16S rRNA)"; 
                                         Is_taxa = true; 
                                         Is_tree = true; 
                                         Is_cp = false;
                                         Is_func = false; 
                                         break;
                                         
                               case 'b': 
                               case 'B':
                               default: DB_id = 'B';
                                        DB_domain = 1;
                                        DB_name = "gg_13"; 
                                        DB_description = "Bacteria (16S rRNA)"; 
                                        Is_taxa = true; 
                                        Is_tree = true; 
                                        Is_cp = true;
                                        Is_func = true; 
                                        break;
                               }
                        if (Is_taxa) Make_taxa();
                        if (Is_tree) Make_tree();
                        if (Is_func) Make_func();
                                                
                  }
             char Get_Id(){
                    return DB_id;
                    }
             
             int Get_Domain(){
                  return DB_domain;
                  }
                    
             string Get_Path(){
                    return DB_path;
                    }                         
             string Get_Description(){
                    return DB_description;
                    }
             //taxa
             int Load_Copy_Number(hash_map <string, float, std_string_hash> & cp_number){
    
                ifstream in_cp_number(DB_cp_number.c_str(), ifstream::in);
                if (!in_cp_number){
                       cerr << "Warning: Cannot open copy number table : " << DB_cp_number << ", copy number normalization is disabled" << endl;
                       return 0;
                       } 
    
                int count = 0;
                string buffer;
                getline(in_cp_number, buffer);
                while(getline(in_cp_number, buffer)){
                                
                                stringstream strin(buffer);
                                string id;
                                float cp_no;
                                strin >> id >> cp_no;
                                
                                if (cp_number.count(id) == 0)
                                                           cp_number[id] = cp_no;
                                else cerr << "Warning, dup id : " << id << endl;
                                
                                count ++;
                                }
    
               in_cp_number.close();
               in_cp_number.clear();
    
               return count;
               }
    
              unsigned int Read_Taxonomy(hash_map<string, string, std_string_hash> & table){
        
                       ifstream infile(DB_taxa_anno.c_str(), ifstream::in);
                       if (!infile){
                          cerr << "Error: Open Taxonomy Annotation file error : " << DB_taxa_anno << endl;
                          exit(0);
                       }
                              
                       unsigned int count = 0;   
                  
                       string buffer;
                       getline(infile, buffer); //title
    
                       while (getline(infile, buffer)){
          
                             if (buffer.size() == 0) continue;
          
                             string id;
                             string taxa;
          
                             int first_table =  buffer.find('\t', 0); //location of first '\t'
                             int last_table = buffer.rfind('\t', buffer.size()-1);//location of last '\t'
          
                             id = buffer.substr(0, first_table);
                             taxa = buffer.substr(last_table + 1, buffer.size() - last_table -1 );
          
                             if (table.count(id) > 0) {
                               
                                                 cerr << "Error: Loading Taxonomy Annotation error : Duplicate ID : " << id << endl;
                                                 exit(0);
                               
                                                 }
                               else (table)[id] = taxa;
          
                               count ++;
          
                               }
    
                     infile.close();
                     infile.clear();
    
                     return count;    
                     }
              //tree
              string Get_Tree_Id(){
                     return DB_tree_id;
                     }
                            
              string Get_Tree_Order(){
                     return DB_tree_order;
                     }
              
              //func
              string Get_Func_Id(){
                     return DB_func_id;
                     }
              string Get_Func(){
                     return DB_func;
                     }
              string Get_Func_Des(){
                     return DB_func_des;
                     }
              string Get_Func_Pw(){
                     return DB_func_pw;
                     }
              string Get_NSTI(){
                     return DB_nsti;
                     }       
              
              bool Get_Is_Cp(){
                   return Is_cp;
                   }
                   
              bool Get_Is_Func(){
                   return Is_func;
                   }
              
      private:
              char DB_id;
              int DB_domain; //0: bacteria; 1: euk
              string DB_name;
              string DB_description;
              string DB_path;
              string DB_cp_number;
              string DB_taxa_anno;
              
              string DB_tree_id;
              string DB_tree_order;
              
              string DB_func_id;
              string DB_func;
              string DB_func_des;
              string DB_func_pw;
              string DB_nsti;
              
              bool Is_taxa;
              bool Is_tree;
              bool Is_cp;
              bool Is_func;
              
              void Make_taxa(){
                   DB_path += "/databases/";
                   DB_path += DB_name;
                   DB_path += "/";
                   DB_cp_number = DB_path + "copy_number.txt";
                   DB_taxa_anno = DB_path + "taxonomy_annotation.txt";                               
                   }
                   
              void Make_tree(){
                   DB_tree_id = DB_path + "tree/id.txt";
                   DB_tree_order = DB_path + "tree/order.txt";
                   };
              void Make_func(){
                   DB_func_id = DB_path + "KO/ko_id.tab";
                   DB_func = DB_path + "KO/ko.tab";
                   DB_func_des = DB_path + "KO/ko_des.tab";
                   DB_func_pw = DB_path + "KO/ko_pw.tab";
                   DB_nsti = DB_path + "otu_nsti.tab";  
                   }
      };


#endif
