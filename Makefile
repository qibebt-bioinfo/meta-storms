CC:=g++
ifneq (,$(findstring Darwin,$(shell uname)))
	exist = $(shell if [ -e '/usr/local/bin/g++-9' ]; then echo "exist"; else echo "notexist"; fi;)
	ifeq ($(exist),exist)
		CC:=g++-9
	else
		CC:=g++-8
	endif
endif
OMP=-fopenmp
HASHFLG=-Wno-deprecated
BUILDFLG=-w -ffunction-sections -fdata-sections -fmodulo-sched -msse
MYSQL=`mysql_config --cflags --libs`

EXE_DMK_OTU=bin/MetaDB-make-otu
EXE_DSC_OTU=bin/MetaDB-search-otu
EXE_DMK_FUNC=bin/MetaDB-make-func
EXE_DSC_FUNC=bin/MetaDB-search-func
EXE_DMK_SP=bin/MetaDB-make-sp
EXE_DSC_SP=bin/MetaDB-search-sp
EXE_DMG=bin/MetaDB-merge
EXE_MNS=bin/MetaDB-parse-mns
EXE_MET=bin/MetaDB-parse-meta
EXE_PQO=bin/MetaDB-parse-qiime-otu

all:
	$(CC) -o $(EXE_DMK_OTU) src/MetaDB_make_otu.cpp $(HASHLIB) $(BUILDFLG) $(OMP)
	$(CC) -o $(EXE_DSC_OTU) src/MetaDB_search_otu.cpp $(HASHLIB) $(BUILDFLG) $(OMP)
	$(CC) -o $(EXE_DMK_FUNC) src/MetaDB_make_func.cpp $(HASHLIB) $(BUILDFLG) $(OMP)
	$(CC) -o $(EXE_DSC_FUNC) src/MetaDB_search_func.cpp $(HASHLIB) $(BUILDFLG) $(OMP)
	$(CC) -o $(EXE_DMK_SP) src/MetaDB_make_sp.cpp $(HASHLIB) $(BUILDFLG) $(OMP)
	$(CC) -o $(EXE_DSC_SP) src/MetaDB_search_sp.cpp $(HASHLIB) $(BUILDFLG) $(OMP)
	$(CC) -o $(EXE_DMG) src/MetaDB_merge.cpp $(HASHLIB) $(BUILDFLG)
	$(CC) -o $(EXE_MNS) src/MetaDB_parse_mns.cpp $(HASHLIB) $(BUILDFLG)
	$(CC) -o $(EXE_MET) src/MetaDB_parse_meta.cpp $(HASHLIB) $(BUILDFLG)
	$(CC) -o $(EXE_PQO) src/MetaDB_parse_qiime_otu.cpp $(HASHLIB) $(BUILDFLG)
clean:
	rm -rf bin/*
