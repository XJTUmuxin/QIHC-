#include<bits/stdc++.h>
using namespace std;
class correcting_seq{
public:
    string seq;
    int seq_len;
    vector<vector<pair<char,int>>> base_count;
    correcting_seq()=default;
    correcting_seq(const string& str){
        seq = str;
        seq_len = seq.size();
        base_count = vector<vector<pair<char,int>>>(seq_len,vector<pair<char,int>>
        {make_pair('A',0),make_pair('C',0),make_pair('G',0),make_pair('T',0),make_pair('N',0)});
    }
    correcting_seq(const correcting_seq& c_seq){
        seq = c_seq.seq;
        seq_len = c_seq.seq_len;
        base_count = c_seq.base_count;
    }
    correcting_seq& operator=(const correcting_seq& c_seq){
        seq = c_seq.seq;
        seq_len = c_seq.seq_len;
        base_count = c_seq.base_count;
        return *this;
    }
};