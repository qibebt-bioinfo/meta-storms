// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// Updated at Dec 19, 2019
// Updated by Xiaoquan Su, Yufeng Zhang
// comp.h
// Compatible with PM3 or above (16S)
// Incompatible for PM1, PM2 and 18S
// Added unweight
// Added HDD model

#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <map>
#include <set>
#include <omp.h>
#include <stdio.h>

#include "comp_dynamic.h"
#include "utility.h"
#include "MetaDB_entry.h"
#include "MetaDB_index.h"
#include "MetaDB_comper.h"
#include "table_format.h"
#include "qsort.h"
#include "hash.h"

#ifndef METADB_H
#define METADB_H

#define BUFF 10000000 

class _MetaDB{
      
      public:
             _MetaDB(){
                       Is_cp_correct = false;
                       Mem_mode = 0;
                       Index_level = DEFAULT_INDEX_LEVEL;                                         
                       };
             
             _MetaDB(char db){
                       Is_cp_correct = false;
                       Mem_mode = 0;
                       Index_level = DEFAULT_INDEX_LEVEL;
                       
                       MetaDB_Index = _MetaDB_Index(db);
                       MetaDB_Comper = _MetaDB_Comper(db);                
                       };
                                                                          
             int Get_Sample_Number(){                 
                 return Samples.size();                 
                 }
                                            
             //Make database
             int Make_Database(const char * databasefile);        
                            
             //Load database
             int Load_Database(const char * databasefile);
             
             //Set prefix for HDD mode
             void Set_Sample_Hdd_Prefix(string p){
             	  Mem_mode = 1;
                  Sample_prefix = p;
                  }
                          
      protected:             
              vector <_MetaDB_Entry> Samples;  
              bool Is_cp_correct;
              int Index_level; //3: family
              int Mem_mode; //0: RAM; 1: HDD  
              string Sample_prefix; //prefix for HDD mode
                                                       
             _MetaDB_Index MetaDB_Index;   
             _MetaDB_Comper MetaDB_Comper;
             			 
			 string Get_Sample_Hdd_Path(string id);                                    
      };

//Make database
int _MetaDB::Make_Database(const char * databasefile){
    
                            ofstream outfile(databasefile, ios::out);
                            if (!outfile){
                               cerr << "Error: Cannot open database file: " << databasefile << endl;
                               return 0;
                               }
                               
                           //output config
                            if (Is_cp_correct) outfile << "T";
                            else outfile << "F";
                            outfile << Index_level << endl; 
                           
                            for (int i = 0; i < Samples.size(); i ++){                                                                              
                                       
                                       //ID
                                       outfile << Samples[i].Get_ID() << endl;
                                       //abd
                                       outfile << "Abd:";
                                       for (int j = 0; j < MetaDB_Comper.Get_dim(); j ++)
                                           if (Samples[i].Get_Abd()[j] > 0)
                                              outfile << "\t" << j << "\t" << Samples[i].Get_Abd()[j];
                                       outfile << endl;
                                       //index
                                       outfile << "Index:";
                                       for (int j = 0; j < MetaDB_Index.Get_IndexN(); j ++)
                                           if (Samples[i].Get_Index_Abd()[j] > 0)
                                              outfile << "\t" << j << "\t" << Samples[i].Get_Index_Abd()[j];
                                       outfile << endl;
                                       //Other
                                       outfile << Samples[i].Get_File() << endl;
                                                                
                                }
                            
                            return Samples.size();
    }


int _MetaDB::Load_Database(const char * databasefile){
				
                 FILE * fptr = fopen(databasefile, "r");                 
                 if (fptr == NULL){
                              cerr << "Error: Cannot open database file: " << databasefile << endl;
                              return 0;
                              }
                 //config                                               
                //while(iter < file_buffer.size()){ //load config
                char * file_buffer = new char [BUFF];
                while(fgets(file_buffer, BUFF-1, fptr) != NULL){
                                       string buffer = file_buffer;
                                       if (buffer.size() == 0) continue;
                                       if (buffer[0] == '#') continue;
                                       if (buffer[0] == 'T')
                                                     Is_cp_correct = true;
                                       else Is_cp_correct = false;            
                                       break;
                                       }               
         
				while(fgets(file_buffer, BUFF-1, fptr) != NULL){
                                       if (file_buffer[strlen(file_buffer)-1] == '\n')
												file_buffer[strlen(file_buffer)-1] = '\0';                                       
                                       string id = file_buffer;
                                       fgets(file_buffer, BUFF-1, fptr);
                                       string abd_buffer = file_buffer;
                                       fgets(file_buffer, BUFF-1, fptr);
                                       string index_abd_buffer = file_buffer;
                                       fgets(file_buffer, BUFF-1, fptr);
                                       string file = "";
                                       string temp;
                                       
                                       //abd                                       
                                       float * abd = new float [MetaDB_Comper.Get_dim()];
                                       memset(abd, 0, MetaDB_Comper.Get_dim() * sizeof(float));
                                       if (Mem_mode ==0){ //RAM mode                                                                                 
                                          stringstream strin_abd(abd_buffer);
                                       
                                          int a_abd_pos;
                                          float a_abd;
                                       
                                          strin_abd >> temp;                                                                            
                                          while(strin_abd >> a_abd_pos){
                                                       strin_abd >> a_abd;
                                                       abd[a_abd_pos] = a_abd;
                                                       } 
                                       }
                                       //index_abd
                                       float * index_abd = new float [MetaDB_Index.Get_IndexN()];
                                       memset(index_abd, 0, MetaDB_Index.Get_IndexN() * sizeof(float));
                                       stringstream strin_index(index_abd_buffer);
                                                                              
                                       int a_index_pos;
                         
						               float a_index;
                                       
                                       strin_index >> temp;
                                       while(strin_index >> a_index_pos){
                                                         strin_index >> a_index;
                                                         index_abd[a_index_pos] = a_index;
                                                         }
									                                         
                                       //add sample
                                       _MetaDB_Entry a_sample;
                                       if (Mem_mode == 0){ //RAM mode
                                          a_sample = _MetaDB_Entry(id, abd, MetaDB_Comper.Get_dim(), index_abd, MetaDB_Index.Get_IndexN(), file);
                                          a_sample.Set_OTU_Count(abd, MetaDB_Comper.Get_dim());  
                                      	}
                                       else //HDD mode
                                            a_sample = _MetaDB_Entry(id, NULL, 0, index_abd, MetaDB_Index.Get_IndexN(), file);
                                       
                                       Samples.push_back(a_sample);
									                                                                                    
                                       delete [] abd;
                                       delete [] index_abd;  
									                                                                              
                                }
				 fclose(fptr);
				 
                 return Samples.size();
                 }

string _MetaDB::Get_Sample_Hdd_Path(string id){
	
	return (Sample_prefix + "/" + id + ".dat");
	
} 

   
#endif
