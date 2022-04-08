#include<bits/stdc++.h>
#include"map_result.h"
#include"correcting_seq.h"
#include"judgeHeterozygosity.h"
string map_result_path="./../map_result/";
string raw_L_path = "./../raw_L/";
string corrected_path = "./../corrected/";
unordered_map<char,int> base_map{{'A',0},{'C',1},{'G',2},{'T',3},{'N',4}};

void correct_Lu(vector<bool>& hete,correcting_seq& cor_seq){
    if(hete.size()!=cor_seq.seq_len){
        cout<<"correct_Lu error"<<endl;
        exit(1);
    }
    srand((unsigned)time(NULL));
    int length = hete.size();
    for(int i=0;i<length;++i){
        if(hete[i]){
            if(cor_seq.seq[i]==cor_seq.base_count[i][0].first || cor_seq.seq[i]==cor_seq.base_count[i][1].first){
                continue;
            }
            else{
                int index = rand()%2;
                cor_seq.seq[i]=cor_seq.base_count[i][index].first;
            }
        }
        else{
            cor_seq.seq[i] = cor_seq.base_count[i][0].first;
        }
    }
}

int main(){
    unordered_map<string,correcting_seq> Lu;
    ifstream ifs(raw_L_path+"Lu.fasta");
    string line;
    while(ifs){
        string name;
        string seq;
        getline(ifs,name);
        name = name.substr(1,name.size()-1);
        getline(ifs,seq);
        correcting_seq cor_seq(seq);
        Lu[name] = cor_seq;
    }
    ifs.close();
    ifs.clear();
    ifs.open(map_result_path+"S_Lu.out");
    getline(ifs,line);
    while(ifs){
        getline(ifs,line);
        if(!line.empty()){
            map_result result(line);
            if(result.matchPercentage>=0.8 && result.qStrand==result.tStrand && result.numIns==result.numDel){
                auto iter = Lu.find(result.tName);
                int index = 0;
                for(int i=result.tStart;i<result.tEnd;++i){
                    while(result.qAlignedSeq[index]=='-')index++;
                    if(index>=result.qAlignedSeq.size())break;
                    iter->second.base_count[i][base_map[result.qAlignedSeq[index]]].second++;
                    index++;
                } 
            }
        }
    }
    ofstream ofs(corrected_path+"corrected_Lu.fasta");
    for(auto &c_seq:Lu){
        vector<bool>hete=judge(c_seq.second,2);
        correct_Lu(hete,c_seq.second);
        ofs<<">"<<c_seq.first<<endl;
        ofs<<c_seq.second.seq<<endl;
    }
    ofs.close();
    ofs.clear();
    return 0;
}