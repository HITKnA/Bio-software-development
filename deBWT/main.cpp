#include<iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <stdlib.h>
#include<math.h>
#include <queue>
#include <bitset>
#include <stdint.h>
#include "x.h"

using namespace std;
using std::bitset;




int main()
{
    ini_ATCG();
    strRef = inputRef("E.coli.fa");
    cout<< "input e.coli.fa"<<endl;

    Kmer2ListLen = inputKmer("mer_counts_dumps.fa");
    cout<<"input kmer counting file"<<endl;
    sort(Kmer2List, Kmer2List+Kmer2ListLen, cmp);
    cout<<"sort kmerList"<<endl;
//check sorted Kmer2List[]
//    ofstream fout;
//    fout.open("Kmer2List");
//    string k;
//    int a;
//    for (int i=0;i<kmerLen;i++)
//    {
//        k = BinTostring(Kmer2List[i],kmer2);
//        a = getKmerNum(Kmer2List[i],kmer2);
//        fout<<k<<" "<<a<<endl;
//    }
//    fout.close();

    KmerIOListSize = checkIO();
    cout<<"buid kmerIOList"<<endl;
//    ofstream fout;
//    fout.open("KmerIOList");
//    int a;
//    string k;
//    for (int i=0;i<IOlen;i++)
//    {
//        k = BinTostring(KmerList[i]<<(64-kmer2*2+2),kmer2);
//        fout<< k<<" ";
//        fout << (KmerIOList[i]>>63)<<" "<<(KmerIOList[i]<<1>>33)<<endl;
//
//    }
//    fout.close();

    buildBCN() ;
    cout<<"build BCN"<<endl;
//    cout << (bc.substr(0,20));
    buildBwt();
    cout << "build BWT"<<endl;
    ofstream fout;
    fout.open("result");
    char cc[6]={'A', 'C', 'G', 'T', '$', 'X'};
    int refLen = strRef.size();
    for (int i=0; i<=refLen; ++i)
    {
        fout<< cc[BWT[i]];
        if ((i+1)%80==0) fout << "\n";
    }
    fout.close();










    return 0;
}
