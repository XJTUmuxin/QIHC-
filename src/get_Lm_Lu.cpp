#include<bits/stdc++.h>
#include"map_result.h"
using namespace std;
string map_result_path = "./../map_result/";
string raw_L_path = "./../raw_L/";
int main(){
    unordered_set<string> Lm;
    unordered_set<string> Lu;
    ifstream ifs(map_result_path+"L_ref.out");
    string line;
    getline(ifs,line);
    while(ifs){
        getline(ifs,line);
        if(!line.empty()){
            map_result result(line);
            string real_name;
            int index = 0;
            for(;index<result.qName.size();++index){
                if(result.qName[index]=='/')break;
            }
            real_name = result.qName.substr(0,index);
            if(result.matchPercentage>=0.87){
                Lm.insert(real_name);
            }
            else{
                Lu.insert(real_name);
            }
        }
    }
    ifs.close();
    ifs.clear();
    ifs.open(raw_L_path+"raw_L.fasta");
    ofstream ofs1(raw_L_path+"Lm.fasta");
    ofstream ofs2(raw_L_path+"Lu.fasta");
    string name;
    string seq;
    while(ifs){
        getline(ifs,line);
        if(line[0]=='>'){
            if(!name.empty()){
                if(Lm.count(name)){
                    ofs1<<">"<<name<<endl;
                    ofs1<<seq<<endl;
                }
                else if(Lu.count(name)){
                    ofs2<<">"<<name<<endl;
                    ofs2<<seq<<endl;
                }
                name = line.substr(1,line.size()-1);
                seq = "";
            }
        }
        else{
            seq += line;
        }
    }
    if(!name.empty()){
        if(Lm.count(name)){
            ofs1<<">"<<name<<endl;
            ofs1<<seq<<endl;
        }
        else if(Lu.count(name)){
            ofs2<<">"<<name<<endl;
            ofs2<<seq<<endl;
        }
    }
    ifs.close();
    ifs.clear();
    ofs1.close();
    ofs1.clear();
    ofs2.close();
    ofs2.clear();
}