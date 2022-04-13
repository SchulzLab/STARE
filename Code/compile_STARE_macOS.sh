working_dir=$(cd "$(dirname "$0")" && pwd -P)

echo "Compiling formatting scripts"
g++ "$working_dir"/ReplaceInvalidChars.cpp -std=c++11 -O3 -o "$working_dir"/ReplaceInvalidChars
g++ "$working_dir"/PSCM_to_PSEM.cpp -std=c++11 -o "$working_dir"/PSCM_to_PSEM
echo "Compiling TRAP"
g++ "$working_dir"/TRAPmulti.cpp -Xpreprocessor -fopenmp -lomp -std=c++11 -O3 -o "$working_dir"/TRAPmulti
echo "Compiling ABC scoring"
g++ "$working_dir"/STARE_ABCpp.cpp "$working_dir"/STARE_MiscFunctions.cpp -Xpreprocessor -fopenmp -lomp -std=c++11 -O3 -o "$working_dir"/STARE_ABCpp
echo "Compiling Gene-TF mapper"
g++ "$working_dir"/TF_Gene_Scorer_Reshape.cpp "$working_dir"/STARE_MiscFunctions.cpp -Xpreprocessor -fopenmp -lomp -std=c++11 -O3 -o "$working_dir"/TF_Gene_Scorer_Reshape
echo "Done!"