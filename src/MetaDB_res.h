// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// Updated at Dec 11, 2018
// Updated by Xiaoquan Su

#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <map>

#include "MetaDB_entry.h"

#ifndef METADB_RES_H
#define METADB_RES_H

class _MetaDB_Res{
      public:
             friend class _MetaDB_search;
             friend class _MetaDB_search_func;
             _MetaDB_Res(){
                       Sim_value = 0;
                       Sample = NULL;
                       Q_OTU_count = 0;
                       }
                       
             _MetaDB_Res(_MetaDB_Entry * s, int q_otu_count, int res_otu_count, float v){
                                   Sim_value = v;
                                   Q_OTU_count = q_otu_count;
                                   Sample = s;                                   
                                   }
             float Get_Res_Sim(){
                    return Sim_value;
                    }
             string Get_Res_Sample_ID(){
                    return Sample->Get_ID();
                    }
             
             _MetaDB_Entry * Get_Sample(){
                           return Sample;
                           }
             int Get_Query_OTU_Count(){
                 return Q_OTU_count;
                 }
             int Get_Res_OTU_Count(){
                 return Res_OTU_count;
                 }
    
            static int Qsort(_MetaDB_Res * res, int begin, int end);
    
      private:
              float Sim_value;
              _MetaDB_Entry * Sample;
              int Res_OTU_count;
              int Q_OTU_count;
             static int Qsort_partition(_MetaDB_Res * res, int begin, int end);
      };

int _MetaDB_Res::Qsort(_MetaDB_Res * res, int begin, int end){
    
    if (begin < end - 1){
        int q = Qsort_partition(res, begin, end);
        Qsort(res, begin, q);
        Qsort(res, q + 1, end);
        }
    return 0;
    }

int _MetaDB_Res::Qsort_partition(_MetaDB_Res * res, int begin, int end){
    
    float x = res[end -1].Sim_value;
    int i = begin - 1;
    for (int j = begin; j < end - 1; j ++)
        if (res[j].Sim_value >= x){
            i ++;
            _MetaDB_Res res_tmp = res[i];
            res[i] = res[j];
            res[j] = res_tmp;
        }
    
    _MetaDB_Res res_tmp = res[i + 1];
    res[i + 1] = res[end - 1];
    res[end - 1] = res_tmp;
    
    return i + 1;
    }

#endif
