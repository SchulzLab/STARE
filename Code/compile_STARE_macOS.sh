g++ ReplaceInvalidChars.cpp -std=c++11 -O3 -o ReplaceInvalidChars
g++ TRAPmulti.cpp -Xpreprocessor -fopenmp -lomp -O3 -o TRAPmulti -std=c++11 -o TRAPmulti
g++ STARE_ABCpp.cpp STARE_MiscFunctions.cpp -Xpreprocessor -fopenmp -lomp -std=c++11 -O3 -o STARE_ABCpp
g++ TF_Gene_Scorer.cpp STARE_MiscFunctions.cpp -Xpreprocessor -fopenmp -lomp -std=c++11 -O3 -o TF_Gene_Scorer
