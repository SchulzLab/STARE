g++ ReplaceInvalidChars.cpp -std=c++11 -o ReplaceInvalidChars
g++ TRAPmulti.cpp -Xpreprocessor -fopenmp -O3 -o TRAPmulti -lomp -std=c++11
g++ STARE_ABCpp_absolute.cpp STARE_MiscFunctions.cpp -Xpreprocessor -fopenmp -std=c++11 -I /usr/local/boost_1_75_0 -o STARE_ABCpp_absolute -lomp
g++ TF_Gene_Scorer_parallel.cpp STARE_MiscFunctions.cpp -Xpreprocessor -fopenmp -std=c++11 -o TF_Gene_Scorer_parallel -lomp
g++ Reshape_toCellTF_perGene.cpp STARE_MiscFunctions.cpp -std=c++11 -o Reshape_toCellTF_perGene
