// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// Updated at June 9, 2018
// Updated by Xiaoquan Su
// MetaDB_make

#include "MetaDB.h"
#ifndef METADB_MAKE_H
#define METADB_MAKE_H

class _MetaDB_make : public _MetaDB{
      
      public:
              //Make database
             int Make_Database(const char * databasefile, string listfile, string listprefix, bool is_cp, int level);
             int Make_Database(const char * databasefile, _Table_Format * tablefile, bool is_cp, int level);
             int Make_Database(const char * databasefile, const char * mdb, int level);    
             int Make_Database_Hdd(const char * mdb); //only make HDD dat
	  
	  private:
	  		int Make_Hdd_File(string id, float * abd, int n);         
      };

//Make database
int _MetaDB_make::Make_Database(const char * databasefile, string listfile, string listprefix, bool is_cp, int level){
                            Is_cp_correct = is_cp;    
                            Index_level = level;                        
                            MetaDB_Index = _MetaDB_Index(Index_level);
                            
                            ofstream outfile(databasefile, ios::out);
                            if (!outfile){
                               cerr << "Error: Cannot open database file: " << databasefile << endl;
                               return 0;
                               }
                            
                            vector <string> ids;
                            vector <string> files;
                            Load_List(listfile.c_str(), files, ids, listprefix);
                                                                         
                            hash_set <string, std_string_hash> id_dup_check;
                            
                            //output config
                            if (Is_cp_correct) outfile << "T";
                            else outfile << "F";
                            outfile << Index_level << endl;
                 
                            for (int i = 0; i < files.size(); i ++){
                                if (id_dup_check.count(ids[i]) != 0){
                                    cerr << "Warning: Duplicated sample ID: " << ids[i] << ", skipped" << endl;
                                    continue;
                                    }
                                //abd
                                float * abd =  new float [Comp_tree.Get_LeafN()];
                                memset(abd, 0, Comp_tree.Get_LeafN() * sizeof(float));
                                Comp_tree.Load_abd(files[i].c_str(), abd, is_cp);
                                //index
                                float * index_abd = new float [MetaDB_Index.Get_IndexN()];
                                memset(index_abd, 0, MetaDB_Index.Get_IndexN() * sizeof(float));
                                MetaDB_Index.Make_Index(abd, index_abd);
                                                                                                
                                //output
                                //ID
                                outfile << ids[i] << endl;
                                //abd
                                outfile << "Abd:";
                                for (int j = 0; j < Comp_tree.Get_LeafN(); j ++)
                                    if (abd[j] > 0)
                                       outfile << "\t" << j << "\t" << abd[j];
                                outfile << endl;
                                if (Mem_mode == 1)
                                	Make_Hdd_File(ids[i], abd, Comp_tree.Get_LeafN());
                                //index
                                outfile << "Index:";
                                for (int j = 0; j < MetaDB_Index.Get_IndexN(); j ++)
                                    if (index_abd[j] > 0)
                                       outfile << "\t" << j << "\t" << index_abd[j];
                                outfile << endl;
                                //Other
                                outfile << ids[i] << endl;

                                id_dup_check.insert(ids[i]);
                                
                                delete [] abd;
                                delete [] index_abd;
                                }
                                
                            outfile.close();
                            outfile.clear();                            
                            return id_dup_check.size();
                            }


int _MetaDB_make::Make_Database(const char * databasefile, _Table_Format * tablefile, bool is_cp, int level){                             
                            Is_cp_correct = is_cp;     
                            Index_level = level;                        
                            MetaDB_Index = _MetaDB_Index(Index_level);
                            
                            ofstream outfile(databasefile, ios::out);
                            if (!outfile){
                               cerr << "Error: Cannot open database file: " << databasefile << endl;
                               return 0;
                               }
                            
                            vector <string> ids = tablefile->Get_Sample_Names();                            
                            hash_set <string, std_string_hash> id_dup_check;
                            
                            //output config
                            if (Is_cp_correct) outfile << "T";
                            else outfile << "F";
                            outfile << Index_level << endl;
                                                        
                            for (int i = 0; i < ids.size(); i ++){                                
                                if (id_dup_check.count(ids[i]) != 0){
                                    cerr << "Warning: Duplicated sample ID: " << ids[i] << ", skipped" << endl;
                                    continue;
                                    }
                                //abd
                                float * abd =  new float [Comp_tree.Get_LeafN()];
                                memset(abd, 0, Comp_tree.Get_LeafN() * sizeof(float));
                                Comp_tree.Load_abd(tablefile, abd, i, is_cp);
                                //index
                                float * index_abd = new float [MetaDB_Index.Get_IndexN()];
                                memset(index_abd, 0, MetaDB_Index.Get_IndexN() * sizeof(float));
                                MetaDB_Index.Make_Index(abd, index_abd);
                                
                                //output
                                //ID
                                outfile << ids[i] << endl;
                                //abd
                                outfile << "Abd:";
                                for (int j = 0; j < Comp_tree.Get_LeafN(); j ++)
                                    if (abd[j] > 0)
                                       outfile << "\t" << j << "\t" << abd[j];
                                outfile << endl;
                                if (Mem_mode == 1)
                                	Make_Hdd_File(ids[i], abd, Comp_tree.Get_LeafN());
                                //index
                                outfile << "Index:";
                                for (int j = 0; j < MetaDB_Index.Get_IndexN(); j ++)
                                    if (index_abd[j] > 0)
                                       outfile << "\t" << j << "\t" << index_abd[j];
                                outfile << endl;
                                //Other
                                outfile << ids[i] << endl;
                                
                                id_dup_check.insert(ids[i]);
                                }   
                                
                            outfile.close();
                            outfile.clear();                            
                            return id_dup_check.size();
                            }    

int _MetaDB_make::Make_Database(const char * databasefile, const char * mdb, int level){   
    
                 Index_level = level;
    
                 ifstream infile(mdb, ios::in);
                 if (!infile){
                              cerr << "Error: Cannot open mdb file: " << mdb << endl;
                              return 0;
                              }
                 
                 ofstream outfile(databasefile, ios::out);
                            if (!outfile){
                               cerr << "Error: Cannot open database file: " << databasefile << endl;
                               return 0;
                               }
                 
                 hash_set <string, std_string_hash> id_dup_check;
                 
                 string buffer;
                 while(getline(infile, buffer)){ //load config
                                       if (buffer.size() == 0) continue;
                                       if (buffer[0] == '#') continue;
                                       if (buffer[0] == 'T')
                                                     Is_cp_correct = true;
                                       else Is_cp_correct = false;                                                                      
                                       break;
                                       }
                                  
                 MetaDB_Index = _MetaDB_Index(Index_level);
                 
                 //output config
                            if (Is_cp_correct) outfile << "T";
                            else outfile << "F";
                            outfile << Index_level << endl;
                                       
                 while(getline(infile, buffer)){//load samples
                           
                                       if (buffer.size() == 0) continue;
                                       if (buffer[0] == '#') continue;
                                       string id;
                                       string file;
                                       string abd_buffer;
                                       string index_abd_buffer;
                                       string temp;
                                       id = buffer;
                                       //abd
                                       getline(infile, abd_buffer);                                       
                                       float * abd = new float [Comp_tree.Get_LeafN()];
                                       memset(abd, 0, Comp_tree.Get_LeafN() * sizeof(float));
                                       stringstream strin_abd(abd_buffer);
                                       
                                       int a_abd_pos;
                                       float a_abd;
                                       
                                       strin_abd >> temp;                                                                            
                                       while(strin_abd >> a_abd_pos){
                                                       strin_abd >> a_abd;
                                                       abd[a_abd_pos] = a_abd;
                                                       }
                     					
                                       //index_abd
                                       getline(infile, index_abd_buffer);
                                       float * index_abd = new float [MetaDB_Index.Get_IndexN()];
                                       memset(index_abd, 0, MetaDB_Index.Get_IndexN() * sizeof(float));
                                       MetaDB_Index.Make_Index(abd, index_abd);
                                                                                                                     
                                       //file path
                                       getline(infile, file);
                                       
                                       if (id_dup_check.count(id) != 0){
                                          cerr << "Warning: Duplicated sample ID: " << id << ", skipped" << endl;
                                          continue;
                                          }
                                       
                                       //output
                                       //ID
                                       outfile << id << endl;
                                       //abd
                                       outfile << "Abd:";
                                       for (int j = 0; j < Comp_tree.Get_LeafN(); j ++)
                                           if (abd[j] > 0)
                                              outfile << "\t" << j << "\t" << abd[j];
                                       outfile << endl;
                                       if (Mem_mode == 1)
                                			Make_Hdd_File(id, abd, Comp_tree.Get_LeafN());
                                       //index
                                       outfile << "Index:";
                                       for (int j = 0; j < MetaDB_Index.Get_IndexN(); j ++)
                                           if (index_abd[j] > 0)
                                              outfile << "\t" << j << "\t" << index_abd[j];
                                       outfile << endl;
                                       //Other
                                       outfile << id << endl;
                                        
                                       id_dup_check.insert(id);
                                        
                                        delete [] abd;
                                        delete [] index_abd;
                                       }
                 
                 infile.close();
                 infile.clear();
                                       
                 outfile.close();
                 outfile.clear(); 
                                       
                 return id_dup_check.size();
                 }       

int _MetaDB_make::Make_Database_Hdd(const char * mdb){  
				
				if (Mem_mode == 0) return 0; 

                 ifstream infile(mdb, ios::in);
                 if (!infile){
                              cerr << "Error: Cannot open mdb file: " << mdb << endl;
                              return 0;
                              }
                 
                 hash_set <string, std_string_hash> id_dup_check;
                 
                 string buffer;
                 while(getline(infile, buffer)){ //load config
                                       if (buffer.size() == 0) continue;
                                       if (buffer[0] == '#') continue;
                                       if (buffer[0] == 'T')
                                                     Is_cp_correct = true;
                                       else Is_cp_correct = false;                                                                      
                                       break;
                                       }
                                       
                 while(getline(infile, buffer)){//load samples
                           
                                       if (buffer.size() == 0) continue;
                                       if (buffer[0] == '#') continue;
                                       string id;
                                       string file;
                                       string abd_buffer;
                                       string index_abd_buffer;
                                       string temp;
                                       id = buffer;
                                       //abd
                                       getline(infile, abd_buffer);                                       
                                       float * abd = new float [Comp_tree.Get_LeafN()];
                                       memset(abd, 0, Comp_tree.Get_LeafN() * sizeof(float));
                                       stringstream strin_abd(abd_buffer);
                                       
                                       int a_abd_pos;
                                       float a_abd;
                                       
                                       strin_abd >> temp;                                                                            
                                       while(strin_abd >> a_abd_pos){
                                                       strin_abd >> a_abd;
                                                       abd[a_abd_pos] = a_abd;
                                                       }
                                                       
                     				   Make_Hdd_File(id, abd, Comp_tree.Get_LeafN());	
                                       //index_abd
                                       getline(infile, index_abd_buffer);
                                                                                                                     
                                       //file path
                                       getline(infile, file);
                                       
                                       if (id_dup_check.count(id) != 0){
                                          cerr << "Warning: Duplicated sample ID: " << id << ", skipped" << endl;
                                          continue;
                                          }

                                        
                                       id_dup_check.insert(id);
                                        
                                       delete [] abd;
                                       }
                 
                 infile.close();
                 infile.clear();
                                       
                 return id_dup_check.size();
                 }  
				                                                                               
int _MetaDB_make::Make_Hdd_File(string id, float * abd, int n){
	
	if (Mem_mode == 0) return 0;
	
	string hdd_path = Get_Sample_Hdd_Path(id);
	ofstream outfile(hdd_path.c_str(), ios::out);
	if (!outfile){
		cerr << "Error: Cannot open file: " << hdd_path << endl;
		return 0;
		}
	int count = 0;
	for (int i = 0; i < n; i ++){
		if (abd[i] > 0){
			outfile  << i << "\t" << abd[i] << endl;
			count ++;
		}
	}
	outfile.close();
	outfile.clear();
	return count;	
}      
#endif
