#include<bits/stdc++.h>
#include"correcting_seq.h"
#include"map_result.h"
#include"judgeHeterozygosity.h"
using namespace std;
string corrected_path = "./../corrected/";
string ref_path = "./../ref/";
string map_result_path = "./../map_result/";
string raw_L_path = "./../raw_L/";
unordered_map<char,int> base_map{{'A',0},{'C',1},{'G',2},{'T',3},{'N',4}};

void correct_Lm(vector<bool>& hete,vector<vector<pair<char,int>>>& base_count,map_result& mr,unordered_map<string,string>& Lm){
    if(hete.size()!=base_count.size()){
        cout<<"correct_Lm error"<<endl;
        exit(1);
    }
    string real_name;
    int index=0;
    for(;index<mr.qName.size();++index){
        if(mr.qName[index]=='/')break;
    }
    real_name = mr.qName.substr(0,index);
    auto iter = Lm.find(real_name);
    srand((unsigned)time(NULL));
    index = mr.tStart;
    for(int i=mr.qStart;i<mr.qEnd;++i){
        if(hete[index]){
            if(iter->second[i]!=base_count[index][0].first || iter->second[i]!=base_count[index][1].first){
                continue;
            }
            else{
                int j = rand()%2;
                iter->second[i] = base_count[index][j].first;
            }
        }
        else{
            iter->second[i] = base_count[index][0].first;
        }
        ++index;
    }
}

int main(){
    //纠错Lm
    string line;
    ifstream ifs(ref_path+"ref.fasta");
    getline(ifs,line);
    getline(ifs,line);
    correcting_seq ref1(line);
    correcting_seq ref2(line);
    ifs.close();
    ifs.clear();  
    ifs.open(map_result_path+"S_ref.out");
    getline(ifs,line);
    while(ifs){
        getline(ifs,line);
        if(!line.empty()){
            map_result result(line);
            if(result.matchPercentage>=0.8 && result.qStrand==result.tStrand && result.numIns==result.numDel){
                int index = 0;
                for(int i=result.tStart;i<result.tEnd;++i){
                    while(result.qAlignedSeq[index]=='-')index++;
                    if(index>=result.qAlignedSeq.size())break;
                    ref1.base_count[i][base_map[result.qAlignedSeq[index]]].second++;
                    index++;
                } 
            }
        }
    }
    vector<bool> hete1 = judge(ref1,2);  //通过S获得的杂合判断结果
    vector<map_result> map_results;
    unordered_map<string,string> Lm;
    ifs.open(raw_L_path+"Lm.fasta");
    while(ifs){
        string name;
        string seq;
        getline(ifs,name);
        name = name.substr(1,name.size()-1);
        getline(ifs,seq);
        if(!seq.empty()){
            Lm[name] = seq;
        }
    }
    ifs.close();
    ifs.clear();
    ifs.open(map_result_path+"Lm_ref.out");
    getline(ifs,line);
    while(ifs){
        getline(ifs,line);
        if(!line.empty()){
            map_result result(line);
            if(result.matchPercentage>=0.87 && result.qStrand==result.tStrand && result.numIns==result.numDel){
                map_results.push_back(result);
                int index = 0;
                for(int i=result.tStart;i<result.tEnd;++i){
                    while(result.qAlignedSeq[index]=='-')index++;
                    if(index>=result.qAlignedSeq.size())break;
                    ref2.base_count[i][base_map[result.qAlignedSeq[index]]].second++;
                    index++;
                } 
            }
        }
    }   
    vector<bool> hete2 = judge(ref2,1);  //通过Lm获得的杂合判断结果
    combine_judge(hete1,hete2);
    combine_base_count(ref1.base_count,ref2.base_count);
    for(auto &v:ref1.base_count){
        sort(v.begin(),v.end(),[](pair<char,int>& a,pair<char,int>& b){return a.second>b.second;});
    }
    ofstream ofs(corrected_path+"corrected_Lm.fasta");
    for(auto &map_result:map_results){
        correct_Lm(hete1,ref1.base_count,map_result,Lm);
        string real_name;
        int index=0;
        for(;index<map_result.qName.size();++index){
            if(map_result.qName[index]=='/')break;
        }
        real_name = map_result.qName.substr(0,index);
        ofs<<">"<<real_name<<":\n";
        ofs<<Lm[real_name]<<"\n";
    }
    ofs.close();
    ofs.clear();
}