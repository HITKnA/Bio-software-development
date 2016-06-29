#include<iostream>
#ifndef x_H
#define x_H
#include<cstring>
#include <stdlib.h>
#include<stdint.h>
#include <bitset>
#include<algorithm>
#include <cstdio>
#include <limits.h>
#include <math.h>
#define numKmer 4864385

using namespace std;
using std::bitset;

string strRef;
int kmer = 20;
int kmer2 = kmer+2;
uint64_t Kmer2List[numKmer+10];
uint64_t KmerIOList[numKmer+10];
uint64_t KmerList[numKmer+10];
int BCNsize,KmerIOListSize;
int Kmer2ListLen;
uint64_t *BCN;
uint64_t beginKmerInRef;
uint8_t *BWT;
string s_A64="AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
string bc;

int ATCGtoBin[85]={0};
char BinTOATCG[4];
void ini_ATCG(){
    ATCGtoBin['A']=0;
    ATCGtoBin['C']=1;
    ATCGtoBin['G']=2;
    ATCGtoBin['T']=3;
    BinTOATCG[0]='A';
    BinTOATCG[1]='C';
    BinTOATCG[2]='G';
    BinTOATCG[3]='T';
}

// kmer*2bit in the front of 64bit is for string
uint64_t stringToBin(string s, int kmer)
{

    char a;
    int i;
    uint64_t bin=0,temp;
    for (i=0;i<kmer;i++)
    {
        a = s[i];
        temp =  ATCGtoBin[a];
        bin |= temp;
        bin <<= 2;
    }
    bin <<= (64-kmer*2-2);
    return bin;
}

    string BinTostring(uint64_t bin,int kmer)
{
    string a;
    a ="AAAAAAAAAAAAAAAAAAAAAA";
    uint64_t b;
    uint64_t c=0;
    b = 3;
    b<<=62;
    int i;
    for(i=0;i<kmer;i++)
    {
        a[i] = BinTOATCG[ ((b&bin)>>62) ];
        bin<<=2;
    }
    return a;
}

string inputRef(string filename)
{
    ifstream readData (filename.c_str());
    string read;
    char buffer2[100];

    if(!readData)
    {
        cout << "Unable to open myfile";
        exit(1); // terminate with error
    }
    readData.getline(buffer2,100);
    while(!readData.eof())
    {
        readData.getline(buffer2, 100);
        read+=buffer2;
    }
    readData.close();
    return read;
}

int getKmerNum(uint64_t bin, int kmer)
{
    uint64_t mask;
    mask = 0xFFFFF;//kmer = 22
    mask = bin & mask;
    return mask;
}


void bitout(uint64_t a)
{
    bitset<64> intbit(a);
    cout<< intbit;
}

int inputKmer(string filename)
{
    ifstream readData(filename.c_str());
    string strKmer,t;
    if (!readData)
    {
        cout<<"unable to open kmerfile";
        exit(1);
    }
    char buffer[30], buffer2[5];
    int numStrKmer, count=0;
    uint64_t temp1, temp2;
    while(!readData.eof())
    {
        readData.getline(buffer2,5);
        readData.getline(buffer,30);
        numStrKmer = buffer2[1];
        numStrKmer -= 48;
        strKmer = buffer;
        temp1 = stringToBin(strKmer,kmer2);
        Kmer2List[count] = temp1;
        temp2 = numStrKmer;
        Kmer2List[count++] |=  temp2;
    }
    return count-1;
}

bool cmp(const uint64_t a, const uint64_t b)
{
    // return a<b;
    return ((a<<2) < (b<<2));
}

uint64_t MASKkmer=-1, MASKin=-1, MASKout=-1, MASKnum=-1;
void iniIOMASK()
{
    MASKkmer = MASKkmer << 2 >> 2 >> (64 - 2*kmer2 + 2) << (64 - kmer2*2 + 2);
    MASKin = MASKin >> 62 << 62;
    MASKout = MASKout >> (64 - 2*kmer2) << (64 - 2*kmer2) << (kmer+1)*2 >> (kmer+1)*2;
    MASKnum = MASKnum << (2*kmer2) >> (2*kmer2);
}


int checkIO()
{

    iniIOMASK();
    uint64_t MASK1 = 1;

    bool isIn, isOut;
    int klen = 0;
    BCNsize = 0;
    uint64_t tempMask,tempNum,temp=0, lastKmer=0, nowKmer=0;

    //get 1st kmer
    for (int i=0; i<kmer;++i)
    {
         temp = temp<<2;
         temp |= ATCGtoBin[strRef[i]];
    }
    beginKmerInRef  =  temp << (64-kmer*2)>>2;
//    bitout(beginKmerInRef);
//    string t;
//    t = BinTostring(beginKmerInRef,kmer2);
//    cout<<endl;
//    cout<<t;

    int i,j,k;
    for (i=0;i<numKmer;)
    {
        isIn = false;
        isOut = false;
        tempMask = 0;
        tempNum = 0;
        j = i+1;
        nowKmer = Kmer2List[i]&MASKkmer;

        if( (nowKmer>=beginKmerInRef && beginKmerInRef>=lastKmer) && (nowKmer == beginKmerInRef))
        {
            isIn = true;
            if ((Kmer2List[i]&MASKout)>>(64-kmer2*2) != ATCGtoBin[strRef[kmer]])
                isOut = true;
            ++tempNum;
        }

        lastKmer = nowKmer;

        while (j<numKmer && ((Kmer2List[j]&MASKkmer)==(nowKmer)))
        {
            if ((Kmer2List[j]&MASKin) != (Kmer2List[i]&MASKin)) isIn = true;
            if ((Kmer2List[j]&MASKout) != (Kmer2List[i]&MASKout)) isOut = true;
            j++;
        }

        if (isIn || isOut){
        if (isOut)
            tempMask |= (MASK1<<63);
        if (isIn)
        {
            for (k=i;k<j;++k)
            {
                tempNum += (Kmer2List[k]&MASKnum);
            }
            tempMask |= BCNsize;
            tempMask |= (tempNum<<32);
            BCNsize += tempNum;
        }
        KmerIOList[klen] = tempMask;
        KmerList[klen] = (Kmer2List[i]&MASKkmer)>>(64-kmer2*2+2);
        klen++;
        }
        i = j;

    }

    return klen;
}

uint64_t BCNmaskOut=-1, BCNmaskIn=-1, BCNmaskIndex=-1,BCNmaskK=-1;
void iniBCNmask()
{
    BCNmaskOut = BCNmaskOut >> 63 << 63;
    BCNmaskIn = BCNmaskIn >> 32 << 32 << 1 >> 1;
    BCNmaskIndex = BCNmaskIndex << 32 >> 32;
    BCNmaskK = BCNmaskK << (64-kmer*2) >> (64-kmer*2);
}

int search_k(uint64_t n)
{
    //binary search
    size_t l=0, r=KmerIOListSize-1, m;
    while (l < r)
    {
        m = (l+r)/2;
        if (n<= KmerList[m]) r = m;
        else l = m + 1;
    }
    if (KmerList[r] == n) return r;
    return -1;
}

int brCodeIndex = 1;
 void buildBCN()
 {
    iniBCNmask();
    BCN = new uint64_t[BCNsize];
    for (int i=0;i<BCNsize;i++) BCN[i]=0;
    uint64_t tmp=0,BCNindex, BCNlength, tk, tempI, sum = 0,temp2;

//    string ss;
//    ss = BinTostring((temp<<64-kmer*2),kmer2);
//    cout<<ss;

    bc = strRef;
    brCodeIndex = 1;
    int refLen = strRef.size();
    for (int i=0, index; i<refLen-1; ++i)
    {
        tmp = (tmp << 2) | ATCGtoBin[strRef[i]];
        if (i>=kmer-1)
        {
            temp2 = tmp&BCNmaskK;
            index = search_k(temp2);

            if (index!=-1)
            {
                tempI = KmerIOList[index];

                //if is multip in
                if (tempI&BCNmaskIn)
                {
                    BCNindex = (tempI&BCNmaskIndex);
                    BCNlength = ((tempI&BCNmaskIn) >> 32);
                    tk = BCNindex;
                    //the BC index stemp2ts from 1
                    while (tk<BCNsize && BCN[tk])
                        ++tk;
                    if (tk < (BCNindex + BCNlength))
                    {
                      if (i>=kmer)
                            BCN[tk] = (brCodeIndex << 3)|ATCGtoBin[strRef[i-kmer]];
                        else
                            BCN[tk] = (brCodeIndex << 3)|4;
                    }
                }
                if (tempI&BCNmaskOut)
                {
                    bc[brCodeIndex++] = strRef[i+1];
                }
            }
        }
    }
    bc.erase(bc.begin()+brCodeIndex, bc.end());

 }

bool bwt_cmp(const uint64_t a, const uint64_t b)
{
    uint64_t ta=(a>>3), tb=(b>>3);
    while (ta<brCodeIndex && tb<brCodeIndex && bc[ta] == bc[tb])
    {
        ++ta;
        ++tb;
    }
    if (ta==brCodeIndex) return true;
    if (tb==brCodeIndex) return false;
    return (bc[ta] < bc[tb]);
}

void buildBwt()
{
    uint64_t mask_code = (-1 << 61 >> 61),last_kmer = 0,now_kmer,BCNindex,BCN_len;
    uint8_t tmp_code;
    int refLen=strRef.size(),bwt_index = 0;
    BWT = new uint8_t[refLen];

    uint64_t *last_string, tmp,tmp_num;
    last_string = new uint64_t[kmer+1];
    for (int i=kmer+1; i>=1; --i)
    {
        for (int j=refLen-i; j<refLen; ++j) tmp = (tmp << 2) | ATCGtoBin[strRef[j]];
        last_string[kmer+1-i] = tmp<<(64-(i)*2);
        // cout << last_string[kmer+1-i] << endl;
    }
    sort(last_string, last_string+kmer+1, cmp);
//
//    for (int i=0;i<kmer;i++)
//    {
//        bitout(last_string[i]);
//        cout<<endl;
//    }

    uint64_t mask_c=-1, mask_l=-1, mask_r=-1, mask_n=-1,mask_index=-1;
    mask_c = mask_c << 2 >> 2 >> (64 - 2*kmer2 + 2) << (64 - kmer2*2 + 2);
    mask_l = mask_l >> 62 << 62;
    mask_r = mask_r >> (64 - 2*kmer2) << (64 - 2*kmer2) << (kmer+1)*2 >> (kmer+1)*2;
    mask_n = mask_n << (2*kmer2) >> (2*kmer2);
    mask_index = mask_index << 32 >> 32;
    uint64_t mask_in = -1 >> 32 << 32 << 1 >> 1;
    uint64_t begin_kmer = beginKmerInRef;

    int j, tmp_index=0, l_index=0, index;
    for (int i=0; i<Kmer2ListLen;)
    {
        //find the kmerIOlist if multi in
        bool is_in=false, is_out=false;
        uint64_t tmp_mask=0;
        j = i+1;
        now_kmer = Kmer2List[i]&mask_c;
        while (l_index < kmer+1 && (last_string[l_index]&mask_c) <= now_kmer)
        {
            BWT[bwt_index++] = ((last_string[l_index++]&mask_l)>>62);
            // BWT[bwt_index-1] = 4;
        }
        if (now_kmer >= begin_kmer && begin_kmer >= last_kmer)
        {
            if (now_kmer == begin_kmer)
            {
                is_in = true;
            }
            else
            {
                BWT[bwt_index++] = 4;
            }
            //else if the begin_kmer no appear in Kmer2List, no need to create the io_info of begin_kmer
        }
        last_kmer = now_kmer;
        while (j<Kmer2ListLen && ((Kmer2List[j]&mask_c) == now_kmer))
        {
            if ((Kmer2List[j]&mask_l) != (Kmer2List[i]&mask_l)) is_in = true;
            if ((Kmer2List[j]&mask_r) !=( Kmer2List[i]&mask_r)) is_out = true;
            ++j;
        }
        if (is_in || is_out)
            tmp_mask = KmerIOList[tmp_index++];
        if (is_in)
        {
            BCNindex = (tmp_mask&mask_index);
            BCN_len = ((tmp_mask&mask_in) >> 32);

            sort(BCN+BCNindex, BCN+BCNindex+BCN_len, bwt_cmp);

            for (int k=BCNindex; k<BCNindex+BCN_len; ++k)
            {
                BWT[bwt_index++] = BCN[k]&mask_code;
            }
        }
        else
        {
            tmp_num=0;
            for (int k = i; k<j; ++k)
            {
                tmp_num += (Kmer2List[k]&mask_n);
            }
            tmp_code = ((Kmer2List[i]&mask_l)>>62);
            while (tmp_num--)
            {
                BWT[bwt_index++] = tmp_code;
            }
        }
        i = j;
    }
}





















#endif
