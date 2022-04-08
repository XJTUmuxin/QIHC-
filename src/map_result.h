#include<bits/stdc++.h>
using namespace std;
class map_result{
public:
    string qName;
    int qLength;
    int qStart;
    int qEnd;
    bool qStrand;
    string tName;
    int tLength;
    int tStart;
    int tEnd;
    bool tStrand;
    int score;
    int numMatch;
    int numMismatch;
    int numIns;
    int numDel;
    int mapQV;
    float matchPercentage;
    string qAlignedSeq;
    string matchPattern;
    string tAlignedSeq;
    map_result(const string &str){
        stringstream strstream(str);
        char temp;
        strstream>>qName>>qLength>>qStart>>qEnd;
        strstream>>temp;
        qStrand = temp=='+'?true:false;
        strstream>>tName>>tLength>>tStart>>tEnd;
        strstream>>temp;
        tStrand = temp=='+'?true:false;
        strstream>>score>>numMatch>>numMismatch>>numIns>>numDel>>mapQV
        >>qAlignedSeq>>matchPattern>>tAlignedSeq;
        matchPercentage = ((float)numMatch/(float)(qEnd-qStart));
        //cout<<qName<<qLength<<qStart<<qEnd<<qStrand<<tName<<tLength<<tStart<<tEnd<<tStrand<<score<<numMatch<<numMismatch<<numIns<<numDel<<mapQV<<qAlignedSeq<<matchPattern<<tAlignedSeq<<endl;
    }
};