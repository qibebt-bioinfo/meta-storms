// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// Updated at Dec 29, 2017
// Updated by Xiaoquan Su
// MetaDB_merge

#include <iostream>
#include <fstream>
#include <sstream>

#include "hash.h"
#include "utility.h"

#ifndef METADB_MERGE_H
#define METADB_MERGE_H

using namespace std;
//Merge
int Merge_Database(const char * databasefile1, const char * databasefile2, const char * databasefile){
    
    ifstream infile1(databasefile1, ios::in);
    if (!infile1){
        cerr << "Error: Cannot open input database file 1: " << databasefile1 << endl;
        return 0;
        }
    
    ifstream infile2(databasefile2, ios::in);
    if (!infile2){
        cerr << "Error: Cannot open input database file 2: " << databasefile2 << endl;
        return 0;
        }
    
    string cfg_1;
    string cfg_2;
    
    string buffer1;
    while(getline(infile1, buffer1)){ //load config
        if (buffer1.size() == 0) continue;
        if (buffer1[0] == '#') continue;
        cfg_1 = buffer1;
        break;
    }
    
    string buffer2;
    while(getline(infile2, buffer2)){ //load config
        if (buffer2.size() == 0) continue;
        if (buffer2[0] == '#') continue;
        cfg_2 = buffer2;
        break;
    }
    
    if (cfg_1 != cfg_2){
        cerr << "Error: 2 database files have different configuration" << endl;
        cerr << "DB1: " << cfg_1 << endl;
        cerr << "DB2: " << cfg_2 << endl;
        return 0;
        }
    
    unsigned count = 0;
    hash_set <string, std_string_hash> id_dup_check;
    
    ofstream outfile(databasefile, ios::out);
    if (!outfile){
        cerr << "Error: Cannot open output database file: " << databasefile << endl;
        return 0;
        }
        
    //config
    outfile << cfg_1 << endl;
    
    //db1
    while(getline(infile1, buffer1)){
        if (buffer1.size() == 0) continue;
        if (buffer1[0] == '#') continue;
        
        string id;
        string abd;
        string index_abd;
        string file_path;
        
        id = buffer1;
        getline(infile1, abd);
        getline(infile1, index_abd);
        getline(infile1, file_path);
        
        if (id_dup_check.count(id) != 0){
            cerr << "Warning: Duplicated sample ID: " << id << ", skipped" << endl;
            continue;
            }
        
        outfile << id << endl;
        outfile << abd << endl;
        outfile << index_abd << endl;
        outfile << file_path << endl;

        id_dup_check.insert(id);
        count ++;
    }
    
    //db2
    while(getline(infile2, buffer2)){
        if (buffer2.size() == 0) continue;
        if (buffer2[0] == '#') continue;
        
        string id;
        string abd;
        string index_abd;
        string file_path;
        
        id = buffer2;
        getline(infile2, abd);
        getline(infile2, index_abd);
        getline(infile2, file_path);
        
        if (id_dup_check.count(id) != 0){
            cerr << "Warning: Duplicated sample ID: " << id << ", skipped" << endl;
            continue;
        }
        
        outfile << id << endl;
        outfile << abd << endl;
        outfile << index_abd << endl;
        outfile << file_path << endl;
        
        id_dup_check.insert(id);
        count ++;
    }
    
    infile1.close();
    infile1.clear();
    
    infile2.close();
    infile2.clear();
    
    outfile.close();
    outfile.clear();
    
    return count;
}

#endif
