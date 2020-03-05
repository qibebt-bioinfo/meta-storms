// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// Updated at Dec 11, 2018
// Updated by Xiaoquan Su

#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <map>
#include <set>

#ifndef METADB_ENTRY_H
#define METADB_ENTRY_H

class _MetaDB_Entry{
      public:
             friend class _MetaDB;
             
             _MetaDB_Entry(){
                            Abd = NULL;
                            Index_abd = NULL;
                            OTU_count = 0;
                            }
              
             _MetaDB_Entry(string id, float * abd, int abd_size, float * index_abd, int index_abd_size, string file){
                                  ID = id;
                                  File = file;
                                  OTU_count = 0;
                                  //abd
                                  if (abd_size > 0){
                                     Abd = new float [abd_size];
                                     memcpy(Abd, abd, abd_size * sizeof(float));
                                     }
                                  
                                  //index
                                  if (index_abd_size > 0){
                                     Index_abd = new float [index_abd_size];
                                     memcpy(Index_abd, index_abd, index_abd_size * sizeof(float));      
                                     }                    
                                  }
             /*
             ~_MetaDB_Entry(){
                              delete [] Abd;
                              delete [] Index_abd;
                              }
             */
             
             int Set_OTU_Count(float * abd, int abd_size){
                 //otu count
                 for (int i = 0; i < abd_size; i ++)
                                      OTU_count += (abd[i] > 0) ? 1: 0;
                 return OTU_count;
                 }
             
             int Get_OTU_Count(){
                 return OTU_count;
                 }
             
             string Get_ID(){
                    return ID;
                    }
             
             float * Get_Abd(){
                   return Abd;
                   }
             
             float * Get_Index_Abd(){
                   return Index_abd;
                   }
             string Get_File(){
                   return File;
                   }
             
      private:
              string ID;
              float * Abd;
              float * Index_abd;
              string File;   
              int OTU_count;                        
      };
#endif
