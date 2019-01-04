CC:=g++
ifneq (,$(findstring Darwin,$(shell uname)))
      CC:=g++-8
endif
HASHLIB=-Wno-deprecated
BUILDFLG=-w -ffunction-sections -fdata-sections -fmodulo-sched -msse
OMP=-fopenmp
EXE_DMK=bin/MetaDB-make
EXE_DSC=bin/MetaDB-search
EXE_DMG=bin/MetaDB-merge
EXE_MNS=bin/MetaDB-parse-mns
EXE_MET=bin/MetaDB-parse-meta
EXE_PQO=bin/MetaDB-parse-qiime-otu

all:
	$(CC) -o $(EXE_DMK) src/MetaDB_make.cpp $(HASHLIB) $(BUILDFLG)
	$(CC) -o $(EXE_DSC) src/MetaDB_search.cpp $(HASHLIB) $(BUILDFLG) $(OMP)
	$(CC) -o $(EXE_DMG) src/MetaDB_merge.cpp $(HASHLIB) $(BUILDFLG)
	$(CC) -o $(EXE_MNS) src/MetaDB_parse_mns.cpp $(HASHLIB) $(BUILDFLG)
	$(CC) -o $(EXE_MET) src/MetaDB_parse_meta.cpp $(HASHLIB) $(BUILDFLG)
	$(CC) -o $(EXE_PQO) src/MetaDB_parse_qiime_otu.cpp $(HASHLIB) $(BUILDFLG)	
clean:
	rm -rf bin/*
