#include<bits/stdc++.h>
using namespace std;
string run_path = "./../run/";
string ref_path = "./../ref/";
string get_ref(string file_name){
    ifstream ifs(file_name);
    string ref;
    while(ifs){
        string temp;
        getline(ifs,temp);
        if(temp[0]=='>')continue;
        ref+=temp;
    }
    ifs.close();
    ifs.clear();
    return ref;
}
int main(){
    string ref;
    string file_name = run_path+"test.contigs.fasta";
    ref = get_ref(file_name);
    ofstream ofs(ref_path+"ref.fasta");
    ofs<<">ref"<<endl;
    ofs<<ref<<endl;
    ofs.close();
    ofs.clear();
}