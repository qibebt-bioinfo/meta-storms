// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// Updated at June 9, 2018
// Updated by Xiaoquan Su
// MetaDB_search

#include <sys/stat.h>
#include <omp.h>

#include "MetaDB.h"
#ifndef METADB_SEARCH_H
#define METADB_SEARCH_H

class _MetaDB_Res{
      public:
             friend class _MetaDB_search;
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

class _MetaDB_search : public _MetaDB{
      public:
             //Single sample
             vector <_MetaDB_Res> Index_Search(const char * queryfile, int n, int f, float t, float min_s, bool is_weight, int coren);
             vector <_MetaDB_Res> Search(const char * queryfile, int n, float min_s, bool is_weight, int coren);
    
             //Sample list
             vector <vector <_MetaDB_Res> > Index_Search(vector <string> queryfile, int n, int f, float t, float min_s, bool is_weight, int coren);
             vector <vector <_MetaDB_Res> > Search(vector <string> queryfile, int n, float min_s, bool is_weight, int coren);
    
             //Table format
             vector <vector <_MetaDB_Res> > Index_Search(_Table_Format * tablefile, int n, int f, float t, float min_s, bool is_weight, int coren);
             vector <vector <_MetaDB_Res> > Search(_Table_Format * tablefile, int n, float min_s, bool is_weight, int coren);  
             
             //Output the result in OTU table
             void Output_OTU_Table(const char * outfilename, const char * queryfile, string queryname, vector <_MetaDB_Res> res);                
             //void Output_OTU_Table(const char * outfilename, const char * queryfile, string queryname, vector <_MetaDB_Res> res, const char * listname); 
             
             
      private:
              int Load_Abd_Hdd(string id, float * abd);                          
              vector <_MetaDB_Res> Index_Search(float * q_abd, int n, int f, float t, float min_s, bool is_weight, int coren);                      
              vector <_MetaDB_Res> Search(float * q_abd, int n, float min_s, bool is_weight, int coren);                        

};

//Search releated methods
//Single sample
vector <_MetaDB_Res> _MetaDB_search::Index_Search(const char * queryfile, int n, int f, float t, float min_s, bool is_weight, int coren){ //top n, index =  f * n
                    
                    float * q_abd = new float [Comp_tree.Get_LeafN()];
                    Comp_tree.Load_abd(queryfile, q_abd, Is_cp_correct);
                    vector <_MetaDB_Res> res = Index_Search(q_abd, n, f, t, min_s, is_weight, coren);
                    delete [] q_abd;
                    return res; 
                    } 

vector <_MetaDB_Res> _MetaDB_search::Search(const char * queryfile, int n, float min_s, bool is_weight, int coren){// top n
                    
                    float * q_abd = new float [Comp_tree.Get_LeafN()];
                    Comp_tree.Load_abd(queryfile, q_abd, Is_cp_correct);
                    vector <_MetaDB_Res> res = Search(q_abd, n, min_s,is_weight, coren);
                    delete [] q_abd;
                    return res;
                    }

//Sample list
vector <vector <_MetaDB_Res> > _MetaDB_search::Index_Search(vector <string> queryfile, int n, int f, float t, float min_s, bool is_weight, int coren){
    
    int sample_num = queryfile.size();

    vector <vector <_MetaDB_Res> > res_return;
    vector <_MetaDB_Res> * res = new vector <_MetaDB_Res> [sample_num];
    
    int single_sample_coren = 1;
    if (sample_num * 2 <= coren)
        single_sample_coren = coren / sample_num;
    
    omp_set_num_threads(coren);
    
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < sample_num; i ++){
        res[i] = Index_Search(queryfile[i].c_str(), n, f, t, min_s, is_weight, single_sample_coren);
        }
    
    for (int i = 0; i < sample_num; i ++)
        res_return.push_back(res[i]);
    
    delete [] res;    
    return res_return;
}

vector <vector <_MetaDB_Res> > _MetaDB_search::Search(vector <string> queryfile, int n, float min_s, bool is_weight, int coren){
    
    int sample_num = queryfile.size();
    
    vector <vector <_MetaDB_Res> > res_return;
    vector <_MetaDB_Res> * res = new vector <_MetaDB_Res> [sample_num];
    
    int single_sample_coren = 1;
    if (sample_num* 2 <= coren)
        single_sample_coren = coren / sample_num;
    
    omp_set_num_threads(coren);
    
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < sample_num; i ++){
        res[i] = Search(queryfile[i].c_str(), n, min_s, is_weight, single_sample_coren);
        }
    
    for (int i = 0; i < sample_num; i ++)
        res_return.push_back(res[i]);
        
    delete [] res;
    return res_return;
    }

//Table format
vector <vector <_MetaDB_Res> > _MetaDB_search::Index_Search(_Table_Format * tablefile, int n, int f, float t, float min_s, bool is_weight, int coren){
    
    int sample_num = tablefile->Get_Sample_Size();
       
    vector <vector <_MetaDB_Res> > res_return;
    vector <_MetaDB_Res> * res = new vector <_MetaDB_Res> [sample_num];
    
    int single_sample_coren = 1;
    if (sample_num * 2 <= coren)
        single_sample_coren = coren / sample_num;
    
    omp_set_num_threads(coren);
    
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < sample_num; i ++){
           float * q_abd = new float [Comp_tree.Get_LeafN()];
           Comp_tree.Load_abd(tablefile, q_abd, i, Is_cp_correct);
           res[i] = Index_Search(q_abd, n, f, t, min_s, is_weight, single_sample_coren);
           delete [] q_abd;
           }
    
    for (int i = 0; i < sample_num; i ++)
        res_return.push_back(res[i]);
    
    delete [] res;   
    return res_return;
    }

vector <vector <_MetaDB_Res> > _MetaDB_search::Search(_Table_Format * tablefile, int n, float min_s, bool is_weight, int coren){
    
    int sample_num = tablefile->Get_Sample_Size();
    
    vector <vector <_MetaDB_Res> > res_return;
    vector <_MetaDB_Res> * res = new vector <_MetaDB_Res> [sample_num];
    
    int single_sample_coren = 1;
    if (sample_num * 2 <= coren)
        single_sample_coren = coren / sample_num;

    omp_set_num_threads(coren);
    
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < sample_num; i ++){
           float * q_abd = new float [Comp_tree.Get_LeafN()];
           Comp_tree.Load_abd(tablefile, q_abd, i, Is_cp_correct);
           res[i] = Search(q_abd, n, min_s, is_weight, single_sample_coren);
           delete [] q_abd;
           }
    
    for (int i = 0; i < sample_num; i ++)
        res_return.push_back(res[i]);
    
    delete [] res;   
    return res_return;
    }

void _MetaDB_search::Output_OTU_Table(const char * outfilename, const char * queryfile, string queryname, vector <_MetaDB_Res> res){
                    
                    vector <string> otus;
                    for (int i = 0; i < Comp_tree.Get_LeafN(); i ++)
                        otus.push_back("OTU_" + Comp_tree.Get_Id(i));
                    
                    _Table_Format table(otus);
                    
                    float * q_abd = new float [Comp_tree.Get_LeafN()];
                    Comp_tree.Load_abd(queryfile, q_abd, Is_cp_correct);
                    
                    vector <float> q_abd_vector(q_abd, q_abd + Comp_tree.Get_LeafN());
                    
                    for (int i = 0; i < q_abd_vector.size(); i ++)
                        q_abd_vector[i] /= 100.0;
                        
                    table.Add_Abd(q_abd_vector, queryname);
                    
                    for (int i = 0; i < res.size(); i ++){
                        vector <float> res_abd_vector;
                        if (Mem_mode == 0){ //RAM                        
                           res_abd_vector = vector <float> ((res[i].Get_Sample())->Get_Abd(), (res[i].Get_Sample())->Get_Abd() + Comp_tree.Get_LeafN());
                           }
                        else{ //HDD 
                              float * res_abd = new float [Comp_tree.Get_LeafN()];
                              Load_Abd_Hdd((res[i].Get_Sample())->Get_ID(), res_abd);                   
                              res_abd_vector = vector <float> (res_abd, res_abd + Comp_tree.Get_LeafN());
                             }
                              
                        for (int j = 0; j < res_abd_vector.size(); j ++)
                        res_abd_vector[j] /= 100.0;
                        table.Add_Abd(res_abd_vector, res[i].Get_Res_Sample_ID());
                    }
                    table.Filter_Empty();
                    table.Output_Table(outfilename);
                    
                    delete [] q_abd;
     
     }      

/*
void _MetaDB_search::Output_OTU_Table(const char * outfilename, const char * queryfile, string queryname, vector <_MetaDB_Res> res, const char * listfilename){
                    
                    set <string> hash;                    
                    //load pcoa list
                    vector <int> pcoa_list;
                    ifstream infile(listfilename, ios::in);
                    if (!infile){
                                 cerr << "Error: Cannot open pcoa list file" << endl;
                                 return;
                                 }
                    
                    int sam_id;
                    while(infile >> sam_id)
                                 pcoa_list.push_back(sam_id);
                    
                    infile.close();
                    infile.clear();
                    
                    vector <string> otus;
                    for (int i = 0; i < Comp_tree.Get_LeafN(); i ++)
                        otus.push_back("OTU_" + Comp_tree.Get_Id(i));
                    
                    _Table_Format table(otus);
                    
                    //load query
                    float * q_abd = new float [Comp_tree.Get_LeafN()];
                    Comp_tree.Load_abd(queryfile, q_abd, Is_cp_correct);
                    
                    vector <float> q_abd_vector(q_abd, q_abd + Comp_tree.Get_LeafN());
                    
                    for (int i = 0; i < Comp_tree.Get_LeafN(); i ++)
                        q_abd_vector[i] /= 100.0;
                        
                    table.Add_Abd(q_abd_vector, queryname);
                                                                                                                                       
                    //load res
                    for (int i = 0; i < res.size(); i ++){
                        vector <float> res_abd_vector((res[i].Get_Sample())->Get_Abd(), (res[i].Get_Sample())->Get_Abd() + Comp_tree.Get_LeafN());
                        for (int j = 0; j < Comp_tree.Get_LeafN(); j ++)
                        res_abd_vector[j] /= 100.0;
                        table.Add_Abd(res_abd_vector, res[i].Get_Res_Sample_ID());
                        hash.insert(res[i].Get_Res_Sample_ID());
                    }
                    
                    //load pcoa
                    for (int i = 0; i < pcoa_list.size(); i ++){
                        if (hash.count(Samples[pcoa_list[i]].Get_ID()) != 0) continue;
                        vector <float> list_abd_vector(Samples[pcoa_list[i]].Get_Abd(), Samples[pcoa_list[i]].Get_Abd() + Comp_tree.Get_LeafN());
                        for (int j = 0; j < Comp_tree.Get_LeafN(); j ++)
                        list_abd_vector[j] /= 100.0;
                        table.Add_Abd(list_abd_vector, Samples[pcoa_list[i]].Get_ID());                        
                        }         
                    
                    table.Filter_Empty();
                    table.Output_Table(outfilename);     
                    
                    delete [] q_abd;
                    }
*/
//private
vector <_MetaDB_Res> _MetaDB_search::Index_Search(float * q_abd, int n, int f, float t, float min_s, bool is_weight, int coren){ //top n, index =  f * n
                    if (n <= 0) n = 1;                    
                    if (n > Samples.size()) n = Samples.size();
                    int cands_size = f * n;
                    if (cands_size > Samples.size()) cands_size = Samples.size();
                    
                    //q_abd count
                    int q_abd_count = 0;
                    for (int i = 0; i < Comp_tree.Get_LeafN(); i ++)
                        q_abd_count += (q_abd[i] > 0) ? 1 : 0;
                    
                    
                    float * q_index_abd = new float [MetaDB_Index.Get_IndexN()];
                    MetaDB_Index.Make_Index(q_abd, q_index_abd);
                    
                    //Index calc
                    int * cands = new int [Samples.size()];
                    float * index_res = new float [Samples.size()];
    
                    omp_set_num_threads(coren);

                    #pragma omp parallel for schedule(dynamic, 1)
                    for (int i = 0; i < Samples.size(); i ++){
                        cands[i] = i;
                        index_res[i] = MetaDB_Index.Calc_Index(q_index_abd, Samples[i].Get_Index_Abd(), t);
                    }
    
                    //index sort
                    qsort(index_res, cands, 0, Samples.size(), cands_size);
    
                    //Index calc end

                    //Res calc
                    _MetaDB_Res * res = new _MetaDB_Res[cands_size];
                    
                    omp_set_num_threads(coren);
                   
                    #pragma omp parallel for schedule(dynamic, 1)
                    for (int i = 0; i < cands_size; i ++){
                        float sim = 0;
                        int res_abd_count = 0;
                        float * sample_abd = NULL;

                        if (Mem_mode == 0){ //RAM
                            sample_abd = Samples[cands[i]].Get_Abd();
                            res_abd_count = Samples[cands[i]].Get_OTU_Count();
                       		}
                        else { //HDD 
                             sample_abd =  new float [Comp_tree.Get_LeafN()];
                             res_abd_count = Load_Abd_Hdd(Samples[cands[i]].Get_ID(), sample_abd);
                             }
                        sim = Comp_tree.Calc_sim(q_abd, sample_abd, is_weight);
                        sim = (sim > 1.0) ? 1.0 : sim;
                        res[i] = _MetaDB_Res(&(Samples[cands[i]]), q_abd_count, res_abd_count, sim);     
                        
                        if (Mem_mode == 1) //HDD
                           delete [] sample_abd;                     
                        }
                    
                    _MetaDB_Res::Qsort(res, 0, cands_size);

                    vector <_MetaDB_Res> res_return;
                    for (int i = 0; i < n; i ++)
                        if (res[i].Sim_value >= min_s)
                           res_return.push_back(res[i]);
                    
                    delete [] q_index_abd;
                    delete [] cands;
                    delete [] index_res;
                    delete [] res;
                    
                    return res_return;
                    } 

vector <_MetaDB_Res> _MetaDB_search::Search(float * q_abd, int n, float min_s, bool is_weight, int coren){// top n
                    if (n <= 0) n = 1;
                    if (n > Samples.size()) n = Samples.size();
                    
                    //q_abd count
                    int q_abd_count = 0;
                    for (int i = 0; i < Comp_tree.Get_LeafN(); i ++)
                        q_abd_count += (q_abd[i] > 0) ? 1 : 0;
                    
                    _MetaDB_Res * res = new _MetaDB_Res[Samples.size()];
                    
                    omp_set_num_threads(coren);
                    
                    #pragma omp parallel for schedule(dynamic, 1)
                    for (int i = 0; i < Samples.size(); i ++){
                        float sim = 0;
                        int res_abd_count = 0;
                        float * sample_abd = NULL;

                        if (Mem_mode == 0){ //RAM
                            sample_abd = Samples[i].Get_Abd();
                            res_abd_count = Samples[i].Get_OTU_Count();
                        	}
                           
                        else { //HDD 
                             sample_abd =  new float [Comp_tree.Get_LeafN()];
                             res_abd_count = Load_Abd_Hdd(Samples[i].Get_ID(), sample_abd);
                             }
                        sim = Comp_tree.Calc_sim(q_abd, sample_abd, is_weight);
                        
                        sim = (sim > 1.0) ? 1.0 : sim;
                        res[i] = _MetaDB_Res(&(Samples[i]), q_abd_count, res_abd_count, sim);     
                        
                        if (Mem_mode == 1)//HDD
                           delete [] sample_abd;                  
                        }
                    
                    //sort
                    _MetaDB_Res::Qsort(res, 0, Samples.size());
    
                    vector <_MetaDB_Res> res_return;
                    for (int i = 0; i < n; i ++)
                    if (res[i].Sim_value >= min_s)
                           res_return.push_back(res[i]);                        
                    
                    delete [] res;
                    
                    return res_return;
                    }
                    
int _MetaDB_search:: Load_Abd_Hdd(string id, float * abd){
					if (Mem_mode == 0) return 0;
					memset(abd, 0, Comp_tree.Get_LeafN() * sizeof(float));
	 				string hdd_path = Get_Sample_Hdd_Path(id);
					 //fgets
                    FILE * fptr = fopen(hdd_path.c_str(), "r");
                    if (fptr == NULL){
                       cerr << "Error: Cannot open input file: " << hdd_path << endl;
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
                   int count = 0;
                   for(int i = 0; i < file_buffer.size(); i ++){
                                       int a_id = 0;
                                       float a_abd = 0;                                       
                                       sscanf(file_buffer[i].c_str(), "%d%f", &a_id, &a_abd);									   
                                       abd[a_id] = a_abd;
                                       count ++;
                                       }                              
					return count;
}

#endif
