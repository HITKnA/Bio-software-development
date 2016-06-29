#include<iostream>
#ifndef rhat_H
#define rhat_H
#include<cstring>
#include <stdlib.h>
#include<stdint.h>
#include <bitset>
#include<algorithm>
#include <cstdio>
#include <limits.h>
#define  pointerlistI  4194304
#define pointerlistJ 20
#define windowMax  10000000
#define readMax 1024
#define readlongmax 25000
#define windowLong 2048
#include "ksw.h"
using namespace std;
using std::bitset;
//void test();

int pointerList [pointerlistI][pointerlistJ]={0};//pow(4,11)
int windowList[windowMax][2]={0};//len(gene)/11
uint32_t binRef[493895];//use 01 to discribe the string ref
int readMaxWindow[1525];
//int readHash[pointerlistI][5];//[hash][start position]
long len1;//ref lenth
string str1;//ref
uint64_t readHash[readlongmax]={0};//32bit readhash 32bit:position
int readRef[readlongmax][3];//[read position][ref position][long=kmer+x]
int **V;
uint32_t *dp;
size_t *path;
size_t p_index ;
ofstream rout;

void test()
{
  cout<<"yeak";
}

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


uint32_t stringToBin(string s)
{

    char a;
    int i;
    uint32_t bin=0,temp;
    for (i=0;i<10;i++)
    {
        a = s[i];
        temp =  ATCGtoBin[a];
        bin |= temp;
        bin <<= 2;
    }
    a = s[10];
    temp = ATCGtoBin[a];
    bin |=  temp;
    return bin;
}

    string BinTostring(uint32_t bin)
{
    string a="AAAAAAAAAAA";
    uint32_t b;
    b = 3;
    int i;
    for(i=0;i<11;i++)
    {
        a[10-i] = BinTOATCG[b&bin];
        bin>>=2;
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

struct node
{
    friend bool operator< (node n1, node n2)
    {
        return n1.window > n2.window;
    }
    int window;
    int index;
};

void bitout(uint32_t bin)
{
    bitset<32> intbit(bin);
    cout<< intbit;
}

int bsearch(uint32_t p_index, int readLen)
{
        //binary search
        uint32_t l=0, r=readLen-1, m;
        while (l < r)
        {
            m = (l+r)/2;
            if (p_index <= (readHash[m]>>32)) r = m;
            else l = m + 1;
        }
        if ((readHash[r]>>32) == p_index) return r;
        return -1;
}

//void getreadpointer(string s)
//{
//    uint32_t current,next = 0;
//    for (size_t i=0; i &lt; s.length(); i++)
//    {
//        next = (next &lt;&lt; 2) | get_c[s[i]];
//        current = next & MASK;
//        if (i &gt;= 10)
//        {
//                read_wlist[i-10] = current&lt;&lt;32|i-10;
//        }
//    }
//}

int readWindow(ifstream &readData, int readn)
{

         char readbuffer[readlongmax],readbuffer2[20];
        string read,check;
        readData.getline(readbuffer2, 20);
        if (readData.eof()) return 0;
        check = readbuffer2;
        readData.getline(readbuffer,readlongmax);
        read = readbuffer;
        readData.getline(readbuffer,20);
        readData.getline(readbuffer, readlongmax);
     int i,j,kmer=11;
        int readLen = read.length();
        if (readLen<11)
        {
            cout<<readn<<endl;
            cout<<check;
        }
     // a 2000 long read to a certain window
    int readList[readMax][3]={0};// readi hash , windows start position, windows end position
    queue <node> q[readMax];
    node x;
    string readi;
    uint32_t readbin,selectBin,temp;
    selectBin = stringToBin("ATTTTTTTTTT");
    readbin = stringToBin(read.substr(readLen/2-kmer,kmer));
    for (i=0;i<readMax;i++)
    {
        temp = ATCGtoBin[read[i+readLen/2]];
        readbin &= selectBin;
        readbin <<= 2;
        readbin |=  temp;
        readList[i][0] = readbin;
        if (pointerList[readList[i][0]][pointerlistJ-1]==0)  continue;
        readList[i][1] = pointerList[readbin][pointerlistJ-2];
        readList[i][2] = pointerList[readbin][pointerlistJ-1];
        for (j=readList[i][1];j<=readList[i][2];j++)
        {
            x.index=i;
            x.window=windowList[j][1];
            q[i].push(x);
        }
    }

    priority_queue <node> qn;
    for (i=0;i<readMax;i++)
    {
        if (q[i].empty()) continue;
        qn.push(q[i].front());
        q[i].pop();
    }



    int a, count = 0, findWindow[readMax*18];
    while(!qn.empty())
    {
        findWindow[count]=qn.top().window;
        count++;
        a = qn.top().index;
        qn.pop();
        if(!q[a].empty())
        {
            qn.push(q[a].front());
            q[a].pop();
        }
    }

    int max_count=0,max_win, findWindowSize=count;
    count=1;
    for (i=0;i<findWindowSize;i++)
    {
        if (findWindow[i]<findWindow[i+1])
        {
            if (max_count<count) { max_count=count; max_win=findWindow[i];}
            count =1;
        }
        else
        {
            count++;
        }
    }
    readMaxWindow[readn] = max_win;
    cout<<read[0]<<" " << max_win<<" "<<max_count<<endl;
    return 0;
}

void create_matrix(int node_i)
    {
        V = new int*[node_i];
        for (int i=0; i<node_i; ++i) V[i] = new int[node_i];
        uint32_t iw, ir, jw, jr;
    	int il;

        for (int i=0; i<node_i-1; ++i)
        {
            iw = readRef[i][1]; //readRef[i][1];
            ir = readRef[i][0]; //readRef[i][0];
            il = readRef[i][2]; //.len;
            for (int j=i+1; j<node_i; ++j)
            {
                V[i][j]=0;
                jw = readRef[j][1]; //wd[j].index_of_W;
                jr = readRef[j][0]; //wd[j].index_of_R;
                if (jw >= il + iw && jr >= il + ir && jr <= ir + il + 1024)
                {
                    V[i][j] = 1;
                }
            }
        }
    }

void find_path(int node_i)
{
    	dp = new uint32_t[node_i];
    	memset(dp, 0, sizeof(dp));
    	for (size_t i=0; i<node_i-1; ++i)
    	{
    		for (size_t j=i+1; j<node_i; ++j)
    		{
    			if (V[i][j]) dp[j] = (dp[j] > readRef[j][2]+dp[i] ? dp[j] : readRef[j][2]+dp[i]);
    		}
    	}
    	size_t thenode = node_i - 1;

    	//path from 1 to p_index) is the index of path string
    	size_t p_tmp;
        p_index = 0;
    	path = new size_t[node_i];
    	while (thenode != 0)
    	{
    		p_tmp = thenode;
    		for (size_t i = 0; i< p_tmp; ++i)
    		{
    			if (i==thenode || V[i][thenode]==0) continue;
    			if (dp[thenode] == readRef[thenode][2] + dp[i])
    			{
    				thenode = i;
    				// PX(i);
    				path[p_index++] = i;
    				break;
    			}
    		}
    		if(thenode == p_tmp) break;
    	}
}



 int8_t  mat[25];
int8_t match=2, mism=-5, gape=2, gapo=1;
const char correspondTable[] = "MIDNSHP=X";
uint32_t *cigar;
int n_cigar = 0;
uint8_t  readqry[2048];
uint8_t  refqry[2048];
uint8_t *readqry_;
uint8_t *refqry_;
size_t PointerListLen = 11;
const uint8_t seq_nt4_tablet[256] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 2, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 2, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

inline void transIntoDec(uint8_t *transtr,char *str, int length)
{
    for (int i=0;i<length;++i) {
        transtr[i] = seq_nt4_tablet[str[i]];
    }
}
double do_alignment(char* ref, size_t window_up, size_t window_down, char*  read, size_t read_l, string& thecigar)
{

    int k=0;
    for (int i=0; i<5; ++i)
    {
        for (int j=0; j<5; ++j)
        {
            if (i<4 && j<4) mat[k++] = i == j? match : mism;
            else mat[k++]=0;
        }
    }

    int read_len;
    int ref_len;
    int w_;

    uint32_t countM = 0;
    char     trans_cigar[500];
    int startPosCigar = 0;

    int qlen = 0;
    int tlen = 0;


    double score=0;

    size_t last_w=window_up + 11 - 1, last_r=11-1;
    size_t i, w, r, l;

    if (p_index < 3) return 0;
    if (p_index >= 3)
    {
        i = path[p_index-2];
        w = readRef[i][1];
        r = readRef[i][0];
        l = readRef[i][2];
        ref_len = w-last_w;
        read_len = r-last_r;

        if (ref_len > PointerListLen / 2)
        {
            size_t tmp = ref_len - PointerListLen/2;
            ref_len = PointerListLen / 2;
            last_w += tmp;
        }
        if (read_len > PointerListLen / 2)
        {
            size_t tmp = read_len - PointerListLen/2;
            read_len = PointerListLen / 2;
            last_r += tmp;
        }

        // cout << string(read, last_r - PointerListLen + 1, read_len) << " " << string(ref, last_w - PointerListLen + 1, ref_len) << endl;
        if (w!=last_w && r!=last_r)
        {
            transIntoDec(readqry,read + last_r - PointerListLen + 1,read_len);
            transIntoDec(refqry,ref + last_w - PointerListLen + 1,ref_len);
            readqry_ = readqry;
            refqry_ = refqry;
            // cout << r << " " << last_r << " " << ref_len << " " << read_len << endl;
            w_ = max(ref_len, read_len);
            n_cigar = 0;
            // score += ksw_global(read_len,readqry_,ref_len,refqry_,5,mat,gapo,gape,w_,&n_cigar,&cigar);
            score += ksw_extend_core(read_len, readqry_, ref_len, refqry_, 5, mat, gapo, gape, 40, read_len , &qlen, &tlen, &cigar, &n_cigar) - read_len;
            if (n_cigar) {
                for (int z=0;z<n_cigar-1;++z) {
                    startPosCigar += sprintf(trans_cigar,"%u%c",cigar[z]>>4,correspondTable[cigar[z]&0xf]);
                    //sams[countSam]._cigar[startPosCigar] = correspondTable[cigar[z]&0xf];
                    //++startPosCigar;
                    thecigar.append(trans_cigar);
                }

                if (correspondTable[cigar[n_cigar-1]&0xf] == 'M') {
                    countM = cigar[n_cigar-1] >> 4;
                    //startPosCigar += sprintf(trans_cigar,"%u%c",countM,'M');
                    //sams[countSam].cigar.append(trans_cigar);
                } else {
                    startPosCigar += sprintf(trans_cigar,"%u%c",cigar[n_cigar-1]>>4,correspondTable[cigar[n_cigar-1]&0xf]);
                    thecigar.append(trans_cigar);
                }

            }
            free(cigar);
        }
        score+=(l)*match;
        // cout << "---" << countM << endl;
        countM+=l;
        // cout << "---" << countM << endl;
        last_w = w-PointerListLen+l+1;
        last_r = r-PointerListLen+l+1;
    }

    for (size_t o = p_index-2; o>0; --o)
    {
        if (o == p_index) continue;
        i = path[o-1];

        w = readRef[i][1];
        r = readRef[i][0];
        l = readRef[i][2];



        ref_len = w-last_w-PointerListLen+1;
        read_len = r-last_r-PointerListLen+1;


        // cout << last_r << " " << r << " " << last_w << " " << w << " " << l << endl;
        // cout << "----\n" << string(read, last_r, read_len) << " " << string(ref, last_w, ref_len) << endl;

        if (w>last_w+PointerListLen-1|| r>last_r+PointerListLen-1)
        {
            if (ref_len > 2048 || read_len > 2048)
            {
                return -INFINITY;
                cout << "hehe\t" << ref_len << " " << read_len << endl;
            }

            transIntoDec(readqry,read + last_r,read_len);
            transIntoDec(refqry,ref + last_w, ref_len);
            readqry_ = readqry;
            refqry_ = refqry;
            w_ = max(ref_len, read_len);
            n_cigar = 0;
            score += ksw_global(read_len,readqry_,ref_len,refqry_,5,mat,gapo,gape,w_,&n_cigar,&cigar);
            if (n_cigar - 1)
            {
                if (correspondTable[cigar[0]&0xf] == 'M') {
                    countM += (cigar[0] >> 4);
                    startPosCigar += sprintf(trans_cigar,"%uM",countM);
                    thecigar.append(trans_cigar);
                }
                else
                {
                    startPosCigar += sprintf(trans_cigar,"%uM",countM);
                    thecigar.append(trans_cigar);
                    startPosCigar += sprintf(trans_cigar,"%u%c",cigar[0]>>4,correspondTable[cigar[0]&0xf]);
                    //sams[countSam]._cigar[startPosCigar] = correspondTable[cigar[z]&0xf];
                    //++startPosCigar;
                    thecigar.append(trans_cigar);
                }
                countM = 0;
                for (int z=1;z<n_cigar-1;++z) {
                    startPosCigar += sprintf(trans_cigar,"%u%c",cigar[z]>>4,correspondTable[cigar[z]&0xf]);
                    //sams[countSam]._cigar[startPosCigar] = correspondTable[cigar[z]&0xf];
                    //++startPosCigar;
                    thecigar.append(trans_cigar);
                }

                if (correspondTable[cigar[n_cigar-1]&0xf] == 'M') {
                    countM = cigar[n_cigar-1] >> 4;
                    //startPosCigar += sprintf(trans_cigar,"%u%c",countM,'M');
                    //sams[countSam].cigar.append(trans_cigar);
                } else {
                    startPosCigar += sprintf(trans_cigar,"%u%c",cigar[n_cigar-1]>>4,correspondTable[cigar[n_cigar-1]&0xf]);
                    thecigar.append(trans_cigar);
                }

            }
            else if (n_cigar>0)
            {
                if (correspondTable[cigar[0]&0xf] == 'M')
                    countM += (cigar[0] >> 4);
                else {
                    startPosCigar += sprintf(trans_cigar,"%uM",countM);
                    thecigar.append(trans_cigar);
                    startPosCigar += sprintf(trans_cigar,"%u%c",cigar[0]>>4,correspondTable[cigar[0]&0xf]);
                    //sams[countSam]._cigar[startPosCigar] = correspondTable[cigar[z]&0xf];
                    //++startPosCigar;
                    thecigar.append(trans_cigar);
                    countM = 0;

                }
            }
            free(cigar);
        }

        score+=(l)*match;
        // cout << "---" << countM << endl;
        countM+=l;
        // cout << "---" << countM << endl;
        last_w = w-PointerListLen+l+1;
        last_r = r-PointerListLen+l+1;
        // outt << readRef[i][1] << " " << readRef[i][0] << " " << readRef[i][2] << " " << read.substr(readRef[i][0] + 1 - PointerListLen, readRef[i][2]) << " " << dna_f.substr(readRef[i][1] + 1 - PointerListLen, readRef[i][2])<< endl;
    }

    // printf("****\n");
    if (p_index >= 3)
    {
        w = window_down;
        r = read_l;

        ref_len = w-last_w;
        read_len = r-last_r;

        if (ref_len > PointerListLen/2) ref_len = PointerListLen / 2;
        if (read_len > PointerListLen/2) read_len = PointerListLen / 2;

        // cout << last_r << " " << r << " " << last_w << " " << w << " " << l << endl;
        // cout << "----\n" << string(read, last_r, read_len) << " " << string(ref, last_w, ref_len) << endl;

        if (w>last_w && r>last_r)
        {
            transIntoDec(readqry,read + last_r,read_len);
            transIntoDec(refqry,ref + last_w, ref_len);
            readqry_ = readqry;
            refqry_ = refqry;
            w_ = max(ref_len, read_len);
            n_cigar = 0;
            qlen = 0;
            tlen = 0;
            score += ksw_extend_core(read_len, readqry_, ref_len, refqry_, 5, mat, gapo, gape, 40, read_len , &qlen, &tlen, &cigar, &n_cigar) - read_len;
            // score += ksw_global(read_len,readqry_,ref_len,refqry_,5,mat,gapo,gape,w_,&n_cigar,&cigar);
            if (n_cigar) {
                if (correspondTable[cigar[0]&0xf] == 'M') {
                    countM += (cigar[0] >> 4);
                    startPosCigar += sprintf(trans_cigar,"%uM",countM);
                    thecigar.append(trans_cigar);
                }
                else
                {
                    startPosCigar += sprintf(trans_cigar,"%uM",countM);
                    thecigar.append(trans_cigar);
                    startPosCigar += sprintf(trans_cigar,"%u%c",cigar[0]>>4,correspondTable[cigar[0]&0xf]);
                    //sams[countSam]._cigar[startPosCigar] = correspondTable[cigar[z]&0xf];
                    //++startPosCigar;
                    thecigar.append(trans_cigar);
                }
                countM = 0;
                for (int z=1;z<n_cigar;++z) {
                    startPosCigar += sprintf(trans_cigar,"%u%c",cigar[z]>>4,correspondTable[cigar[z]&0xf]);
                    //sams[countSam]._cigar[startPosCigar] = correspondTable[cigar[z]&0xf];
                    //++startPosCigar;
                    thecigar.append(trans_cigar);
                }

            }
            free(cigar);
        }
    }

    // printf("%lf\n", score);

}

int readWindow_array(string read, int readn, const char* c_ref)
{


    int i,j,kmer=11;
    int readLen = read.length();
     // a 2000 long read to a certain window
    int readList[readMax][3]={0};// readi hash , windows start position, windows end position
    node x;
    string readi;
    uint32_t readbin,selectBin,temp;
    selectBin = stringToBin("ATTTTTTTTTT");
    readbin = stringToBin(read.substr(readLen/2-kmer,kmer));
    int COUNT[4824]={0};
    for (i=0;i<readMax;i++)
    {
        temp = ATCGtoBin[read[i+readLen/2]];
        readbin &= selectBin;
        readbin <<= 2;
        readbin |=  temp;
        readList[i][0] = readbin;
        if (pointerList[readList[i][0]][pointerlistJ-1]==0)  continue;
        readList[i][1] = pointerList[readbin][pointerlistJ-2];
        readList[i][2] = pointerList[readbin][pointerlistJ-1];
        for (j=readList[i][1];j<=readList[i][2];j++)
        {
                COUNT[windowList[j][1]] ++;
            }
    }

    int max_count=0,max_win, count;
    count=0;

    for(i=0;i<4824;i++)
    {
       if (COUNT[i]>max_count)
        {
            max_count = COUNT[i];
            max_win = i;
        }
    }
    readMaxWindow[readn] = max_win;
    //cout<<read[0]<<" " << max_win<<" "<<max_count<<endl;

    memset(readHash,0,sizeof(readHash));
    readbin = stringToBin(read.substr(0,kmer));
    readHash[0] |= readbin;
    readHash[0]<<=32;
    for (i=kmer;i<readLen;i++)
    {
        temp = ATCGtoBin[read[i]];
        readbin &=  selectBin;
        readbin <<= 2;
        readbin |= temp;
        readHash[i-kmer+1] |= readbin;
        readHash[i-kmer+1]<<=32;
        temp = i-kmer+1;
        readHash[i-kmer+1] |= temp;
    }
    sort(readHash, readHash+readLen-kmer+1);

    //get window string start and end
    string refString;
    int windowStart = max_win * (windowLong/2);
    int readLeft = readLen/2 - kmer;
    if (readLeft<0) readLeft=0;
    windowStart -= readLeft;
    if (windowStart<0) windowStart=0;
    int windowEnd = windowStart +readLen+windowLong/2;
    if (windowEnd>len1) windowEnd = len1;
    refString = str1.substr(windowStart,windowEnd-windowStart);
    int refStringLen = windowEnd-windowStart;

    //scan the ref and find the same lmer on the read
    int position,rfcount=0;
    memset(readRef,0,sizeof(readRef));
    readbin = stringToBin(refString.substr(0,kmer));
    position = bsearch(readbin,readLen);
    if (position!=-1){
        readRef[rfcount][0] = readHash[position]&0xFFFFFFFF;
        readRef[rfcount][1] = 0+windowStart;
        readRef[rfcount++][2] = 0;}
   for (i=1;i<refStringLen-kmer+1;i++)
   {
        readbin &= selectBin;
        readbin <<= 2;
        temp = ATCGtoBin[refString[i+kmer-1]];
        readbin |= temp;
        position = bsearch(readbin,readLen);
        if (position==-1)//can not find lmer in read
        { continue;}
        else
        {
             readRef[rfcount][0] = readHash[position]&0xFFFFFFFF;
             readRef[rfcount][1] = i+windowStart;
             if ((readRef[rfcount][0]-readRef[rfcount-1][0]-readRef[rfcount-1][2] == 1 ) && (readRef[rfcount][1]-readRef[rfcount-1][1]-readRef[rfcount-1][2] == 1) )//extend the lmer
             {
                    readRef[rfcount-1][2] ++;
             }
             else
             {
                    readRef[rfcount++][2] = 0;
             }
        }
   }




    //create graph
    create_matrix(rfcount);
    find_path(rfcount);

    string ss;
    //ofstream fout("t_fill_the_gap");
    const char *c_read = read.c_str();
    do_alignment(const_cast<char*>(c_ref), windowStart, windowStart+refStringLen, const_cast<char*>(c_read), readLen, ss);
    rout<<ss<<endl;
    rout<<read<<endl;
    //out.close();

        //output readhash table
    //ofstream fout("read_hashtable_checkposition");
//    string t_string;
//    for (i=0;i<readLen-kmer+1;i++)
//    {
//        temp = readHash[i]>>32;
//         if ((temp==readHash[i-1]>>32)){
//            t_string = BinTostring(temp);
//            temp = readHash[i]&0xFFFFFFFF;
//            cout<<"the no. of read:"<<readn<<" posit"<<temp<<" "<<t_string<<endl;
//            }
        //fout<<temp<<" "<<BinTostring(temp)<<" ";
        //temp = readHash[i]&0xFFFFFFFF;
        //fout<<temp<<endl;
//    }
    //fout.close();

//    position = bsearch(stringToBin("AAAAACTGTTT"),readLen);
//    position = readHash[position]&0xFFFFFFFF;
//    cout<<position;


    return max_win;
}




#endif
