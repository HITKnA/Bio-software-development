#include<iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <stdlib.h>
#include<math.h>
#include <queue>
#include "rhat.h"
#include "ksw.h"
#include <bitset>
#include <stdint.h>
#define  pointerlistI  4194304
#define pointerlistJ 20
#define windowMax  10000000

using namespace std;
using std::bitset;




void outputPointerList(int kmer)
{
    ofstream fout;
    fout.open("output_pointer_bin.txt");
    int i,j;
    for  (i=0;i<pointerlistI;i++)
    {
        if (pointerList[i][pointerlistJ-2]==0) continue;
        fout <<BinTostring(i)<<" ";
       // for (j=0;j<pointerlistJ;j++)
        {
            //if (pointerList[i][j]==0) break;
            //fout<<pointerList[i][j]<<" ";
            j = pointerList[i][pointerlistJ-2];
            for (j;j<=pointerList[i][pointerlistJ-1];j++)
            {
                fout <<windowList[j][1]<<" ";
            }
        }
        fout<<endl;
    }
    fout.close();
}



   int main(int argc, const char * argv[]){


    //input file


    str1 = inputRef("test.fa");
    len1 = str1.length();
    cout << len1<<endl;
    cout << "The origin sequence is :   "<<str1[1]<<endl;


    int i, kmer=11,plIndex,j=0;
    string kstring;
    uint32_t kbin,temp,selectBin;
    bitset<32> intbit;
    kstring = str1.substr(0,kmer);

    ini_ATCG();
    kbin = stringToBin(kstring);
    pointerList[kbin][0] = kmer-1;
    intbit |= kbin;

    selectBin = stringToBin("ATTTTTTTTTT");
    int brcount=0;
    binRef[brcount++] = kbin;
    for (i=kmer;i<len1;i++ )
    {
        temp = ATCGtoBin[str1[i]];
        kbin &=  selectBin;
        kbin <<= 2;
        kbin |= temp;
        while (pointerList[kbin][j] != 0)
        {j+=1;}
        pointerList[kbin][j] = i;
        j = 0;
        if ((i+1)%kmer==0)
            binRef[brcount++]=kbin;//translate string to bin
    }

    int wlIndex=0;

    for (i=0;i<pointerlistI;i++)
    {
        if (pointerList[i][0]==0) continue;
        for (j=0;j<pointerlistJ-2;j++)
        {
            if (pointerList[i][j]==0) break;
            windowList[wlIndex][0] = i;
            windowList[wlIndex][1] = pointerList[i][j]/(windowLong /2);
            wlIndex += 1;
            windowList[wlIndex][0] = i;
            windowList[wlIndex][1] = pointerList[i][j]/(windowLong /2)+1;
            wlIndex += 1;
        }
    }

    pointerList[windowList[0][0]][pointerlistJ-2]=0;//find windowList index in pointerList
    for (i=1;i<windowMax;i++)
    {
         if  (windowList[i][0]>windowList[i-1][0])
         {
            pointerList[windowList[i][0]][pointerlistJ-2] = i;
            pointerList[windowList[i-1][0]][pointerlistJ-1] = i-1;
         }
    }
   //outputPointerList(kmer);

    ifstream readData ("readtest");

    if(!readData)
    {
        cout << "Unable to open myfile";
        exit(1); // terminate with error
    }
    char readbuffer[readlongmax],readbuffer2[20];
    string read;
    rout.open("resultx");
    const char *c_ref = str1.c_str();
    for (i=1;i<1250;i++){
        readData.getline(readbuffer2, 20);
        if (readData.eof()) break;
        readData.getline(readbuffer,readlongmax);
        read = readbuffer;
        readData.getline(readbuffer,20);
        readData.getline(readbuffer, readlongmax);
         j = i;
        rout<<"read"<<j<<endl;
        readWindow_array(read,i,c_ref);
//        readData.getline(readbuffer2, 20);
//        if (readData.eof()) break;
//        readData.getline(readbuffer,readlongmax);
//        read = readbuffer;
//        readData.getline(readbuffer,20);
//        readData.getline(readbuffer, readlongmax);

        //memset(readHash,0,sizeof(readHash));

    }
    readData.close();
    rout.close();








   return 0;
  }
