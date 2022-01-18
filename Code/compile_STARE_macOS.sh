working_dir=$(cd "$(dirname "$0")" && pwd -P)

g++ "$working_dir"/ReplaceInvalidChars.cpp -std=c++11 -O3 -o "$working_dir"/ReplaceInvalidChars
g++ "$working_dir"/TRAPmulti.cpp -Xpreprocessor -fopenmp -lomp -O3 -std=c++11 -o "$working_dir"/TRAPmulti
g++ "$working_dir"/STARE_ABCpp.cpp "$working_dir"/STARE_MiscFunctions.cpp -Xpreprocessor -fopenmp -lomp -std=c++11 -O3 -o "$working_dir"/STARE_ABCpp
g++ "$working_dir"/TF_Gene_Scorer.cpp "$working_dir"/STARE_MiscFunctions.cpp -Xpreprocessor -fopenmp -lomp -std=c++11 -O3 -o "$working_dir"/TF_Gene_Scorer
