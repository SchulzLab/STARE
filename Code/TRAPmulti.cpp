//Helge Roider
//Max Planck Institute for Molecular Genetics - Berlin
//DATE: 16/06/2008

//Adapted by Florian Schmidt and Marcel H. Schulz
//Cluster of Excellence on Multimodal Computing and Interaction and Max Planck Insitute for Informatics - Saarbruecken
//DATE: 17/11/2016

// Adapted by Dennis Hecker August 2021, reduced loops in the affinity calculation, changed sequence reading to stream
// buffer. Institute of Cardiovascular Regeneration - Goethe University Frankfurt

/*
* MacOS:
* g++ TRAPmulti.cpp -Xpreprocessor -fopenmp -lomp -O3 -std=c++11 -o TRAPmulti
* Linux:
* g++ TRAPmulti.cpp -fopenmp -O3 -std=c++11 -o TRAPmulti
*/

#include<stdio.h>
#include<cmath>
#include<stdlib.h>
#include<fstream>
#include<iostream>
#include<string.h>
#include<iomanip>
#include<omp.h>

using namespace std;

const int A = 0;
const int C = 1;
const int G = 2;
const int T = 3;


int main(int argc, char *argv[]){

    if(argc < 4){
        cerr << "\nINPUT:\n\t1. ENERGY MATRIX (from PSCM_to_PSEM)\n\t2. FASTA FILE\n\t3. File with the maximal sequence length\n\t4. NUMBER OF CORES TO BE USED";
        exit(1);
    }
    string transfac_file = argv[1];
    string fasta_file = argv[2];
    string max_sequence_file = argv[3];
    unsigned int num_workers=atoi(argv[4]);
    
    //----------------------------------------------------------------
    //Initialise VARIABLES
    //----------------------------------------------------------------
    string row;
    ifstream seq_length(max_sequence_file);
    int max_seq_length;
    if(seq_length){
        getline(seq_length, row);
        max_seq_length = stoi(row);
    }
    else {
        cerr << "Sequence length file not opened\n" << max_sequence_file;
        exit(1);
    }

    int maxMotifLength = 0;
    int numOfFactors = -1;
    ifstream transfac1(transfac_file);
    int m = 0;
    if(!transfac1){
        cerr << "Matrix file not opened\n" << transfac_file;
        exit(1);
    }

    while(!transfac1.eof()){
        getline(transfac1,row);
        if((row.substr(0,1) == "#")||(row.substr(0,1) == "")){
            continue;
        }
        if(row.substr(0,1) == ">"){
            numOfFactors++;
            if(m > maxMotifLength){
                maxMotifLength = m;
            }
            m = 0;
            continue;
        }
        m++;
    }

    if(m > maxMotifLength){
        maxMotifLength = m;
    }
    transfac1.close();

    maxMotifLength++;
    numOfFactors++;

    double lnR0[numOfFactors];
    string factornames[numOfFactors];
    int motiflength[numOfFactors];

    double *** complement;
    double *** matrix;
    matrix = new double ** [numOfFactors];
    complement = new double ** [numOfFactors];
    for(int i = 0; i < numOfFactors; i++){
        matrix[i] = new double * [maxMotifLength];
        complement[i] = new double * [maxMotifLength];
        for(int j = 0; j < maxMotifLength; j++){
            matrix[i][j] = new double[4];
            complement[i][j] = new double[4];
        }
    }


    //----------------------------------------------------------------
    //READ POSITION WEIGHT MATRICES
    //----------------------------------------------------------------

    //reading file variables
    int factors = -1;
    string word[100]; //elements in row
    double max = 0; //consensus base count
    string delimiters = " \t"; //word seperators in each line
    int reading;
    int start, end;

    ifstream transfac(transfac_file);

    if(!transfac){
        cerr << "Matrix file not opened\n";
        exit(1);
    }

    while(!transfac.eof()){
        getline(transfac,row);
        start = row.find_first_not_of(delimiters);
        if(row.substr(0,1) == "#"){ // commentary line
            continue;
        }
        if(row.substr(0,1) == ""){ // empty line
            continue;
        }
        int i = 0;
        while(start != string::npos){ //split row into tokens - word[]
            end = row.find_first_of(delimiters,start + 1);
            if(end == string::npos){
                end = row.length();
            }
            word[i] = row.substr(start,end - start);
            i++;
            start = row.find_first_not_of(delimiters,end + 1);
        }

        if(i == 0){ // line without content
            continue;
        }

        if(word[0].substr(0,1) == ">"){ //new matrix reached
            factors++;
            motiflength[factors] = 0;
            factornames[factors] = word[1];//word[0].substr(2);//1
            for(int w = 0; w < i; w++){
                if(word[w] == "lnR0:"){
                    lnR0[factors] = strtod(word[w + 1].c_str(),NULL);
                }
            }
            continue;
        }
        matrix[factors][motiflength[factors]][A] = strtod(word[0].c_str(),NULL);
        matrix[factors][motiflength[factors]][C] = strtod(word[1].c_str(),NULL);
        matrix[factors][motiflength[factors]][G] = strtod(word[2].c_str(),NULL);
        matrix[factors][motiflength[factors]][T] = strtod(word[3].c_str(),NULL);
        motiflength[factors]++;
    }
    transfac.close();

    for(int f = 0; f <= factors; f++){
        for(int p = 0; p < motiflength[f]; p++){
            complement[f][motiflength[f]-p-1][A] = matrix[f][p][T];
            complement[f][motiflength[f]-p-1][T] = matrix[f][p][A];
            complement[f][motiflength[f]-p-1][G] = matrix[f][p][C];
            complement[f][motiflength[f]-p-1][C] = matrix[f][p][G];
        }
    }

    //----------------------------------------------------------------
    //READ FASTA FILE
    //----------------------------------------------------------------
    for(int i = 0; i <= factors; i++){
        cout << "\t" << factornames[i];
    }
    cout << "\n";

    string seqID;
    int seqlength;
    string newbases; //sequence in each row of fasta file
    FILE *ReadSequence = fopen(fasta_file.c_str(), "rb");
    char sequence_buffer[max_seq_length];
    if(!ReadSequence){
        cerr << "FASTA file not opened:\n" << fasta_file;
        exit(1);
    }
    while (!feof(ReadSequence)) {
        if (fgets(sequence_buffer, max_seq_length, ReadSequence) != NULL) {
            newbases = sequence_buffer;
            newbases = newbases.substr(0, newbases.size() - 1);  // Removing the trailing \n.
            if (newbases.substr(0, 1) == "#") { //
                continue;
            }
            if (newbases.substr(0, 1) == ">") { //NEW SEQUENCE HEADER
                seqID = newbases.substr(1);
                continue;
            }

            if (newbases.substr(0, 1) != ">" and newbases.size() > 1) { //SEQUENCE LINE
                seqlength = newbases.length();

                cout << seqID;
                //New: initialize array for all affinity values of one TF
                double affinityValuesTF[factors];
                //LOOP OVER FACTORS
#pragma omp parallel for num_threads(num_workers)
                for (int f = 0; f <= factors; f++) {
                    affinityValuesTF[f] = 0;
                    int illegalBase, BASE;
                    double P_combined = 0; //only palindrome correction, for entire seq
                    double product, P_bound_F, P_bound_C;
                    double dE_forward, dE_compl;
                    for (int n = 0; n < seqlength - motiflength[f] + 1; n++) { //LOOP OVER SEQUENCE
                        dE_forward = 0;
                        dE_compl = 0;
                        for (int m = 0; m < motiflength[f]; m++) { //LOOP OVER MOTIF
                            illegalBase = 0;
                            switch (newbases[n + m]) {
                                case 'A':
                                    BASE = 0;
                                    break;
                                case 'C':
                                    BASE = 1;
                                    break;
                                case 'G':
                                    BASE = 2;
                                    break;
                                case 'T':
                                    BASE = 3;
                                    break;
                                case 'a':
                                    BASE = 0;
                                    break;
                                case 'c':
                                    BASE = 1;
                                    break;
                                case 'g':
                                    BASE = 2;
                                    break;
                                case 't':
                                    BASE = 3;
                                    break;
                                default:
                                    illegalBase = 1;
                            }
                            if (illegalBase == 1) { break; }
                            dE_forward += matrix[f][m][BASE];
                            dE_compl += complement[f][m][BASE];
                        }//loop over motif

                        //CALCULATE P(BOUND) FOR CURRENT SITE
                        if (illegalBase == 0) {
                            product = exp(lnR0[f] - dE_forward);
                            P_bound_F = product / (1 + product);
                            product = exp(lnR0[f] - dE_compl);
                            P_bound_C = product / (1 + product);
                            P_combined += P_bound_F + (1 - P_bound_F) * P_bound_C;
                        }
                    }//loop over sequence
                    affinityValuesTF[f] = P_combined;
                }//loop over factors

                //Print the affinity values for the current sequence of all TFs
                for (int f = 0; f <= factors; f++) {
                    cout << "\t" << setprecision(10) << affinityValuesTF[f];
                }
                cout << "\n";
            }
        }
    }//loop over fasta file
    fclose(ReadSequence);

    return 0;
}
