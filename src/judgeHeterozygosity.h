#pragma once
#include<bits/stdc++.h>
#include"map_result.h"
#include"correcting_seq.h"
#define p1 0.2      //第三代测序错误率
#define p2 0.05     //第二代测序错误率
#define p3 0.2
#define p4 0.8
#define p5 0.5
double combine_num(int n,int m){
    if (n == m || m == 0) return 1;

    vector<double> dp(m+1);
    for (int i = 0; i <= n; i++){
        for (int j = min(i, m); j >= 0; j--){
            if (i == j || j == 0) dp[j] = 1;
            else dp[j] = dp[j] + dp[j-1];
        }
    }
    return dp[m];
}

void combine_base_count(vector<vector<pair<char,int>>>& base_count1,vector<vector<pair<char,int>>>& base_count2){
    if(base_count1.size()!=base_count2.size()){
        cout<<"combine_base_count_err"<<endl;
        exit(1);
    }
    for(int i=0;i<base_count1.size();++i){
        sort(base_count1[i].begin(),base_count1[i].end(),[](const pair<char,int>& pci1,const pair<char,int>& pci2){return pci1.first<pci2.first;});
        sort(base_count2[i].begin(),base_count2[i].end(),[](const pair<char,int>& pci1,const pair<char,int>& pci2){return pci1.first<pci2.first;});
        for(int j=0;j<5;++j){
            base_count1[i][j].second+=base_count2[i][j].second;
        }
    }
}

void combine_judge(vector<bool>& judge1,vector<bool>& judge2){
    int len1 = judge1.size();
    int len2 = judge2.size();
    if(len1!=len2){
        exit(1);
    }
    for(int i=0;i<len1;++i){
        if(!judge1[i] || !judge2[i]){
            judge1[i] = false;
        }
    }
}

//mode=1 对应Lm的杂合性判断 mode=2 对应Lu的杂合性判断
vector<bool> judge(correcting_seq& cor_seq,int mode){
    if(mode!=1 && mode!=2){
        cout<<"模式错误"<<endl;
        exit(1);
    }
    vector<bool> judge_result(cor_seq.seq_len);
    vector<int> rd(cor_seq.seq_len,0);              //该位点匹配到的序列数
    vector<int> dl(cor_seq.seq_len,0);              //该位点匹配到的序列对应位点出现的碱基种类
    int index = 0;
    for(auto &v:cor_seq.base_count){
        sort(v.begin(),v.end(),[](pair<char,int>& a,pair<char,int>& b){return a.second>b.second;});
        int temp_rd=0;
        int temp_dl=0;
        for(auto &p:v){
            temp_rd+=p.second;
            if(p.second>0)temp_dl++;
        }
        rd[index] = temp_rd;
        dl[index] = temp_dl;
        index++;
    }
    for(int i=0;i<cor_seq.seq_len;++i){
        double p_err; 
        double p_hete;
        double p_homo;
        int Rd=rd[i];
        if(mode==1){    //Lm
            p_err = p1;
        }
        else{           //Lu
            p_err = p2;
        }
        if(dl[i]<=1){
            judge_result[i] = false;
        }
        else if(dl[i]==2){
            int r1 = cor_seq.base_count[i][0].second;
            int r2 = cor_seq.base_count[i][1].second;
            p_homo = combine_num(Rd,r1)*pow(1-p_err,r1)*pow(p_err,Rd-r1)*combine_num(Rd-r1,r2)*pow(p_err,r2)*pow(1-p_err,Rd-r1-r2)*p3;
            p_hete = combine_num(Rd,r1)*pow(p5,r1)*pow(1-p5,Rd-r1)*combine_num(Rd-r1,r2)*pow(p5,r2)*pow(1-p5,Rd-r1-r2)*p4;

            if(p_homo>p_hete){
                judge_result[i] = false;
            }
            else{
                judge_result[i] = true;
            }
        }
        else if(dl[i]==3){
            int r1 = cor_seq.base_count[i][0].second;
            int r2 = cor_seq.base_count[i][1].second;
            int r3 = cor_seq.base_count[i][2].second;
            p_homo = combine_num(Rd,r1)*pow(1-p_err,r1)*pow(p_err,Rd-r1)*combine_num(Rd-r1,r2)*pow(p_err/2,r2)*pow(1-p_err/2,Rd-r1-r2)*
            combine_num(Rd-r1-r2,r3)*pow(p_err/2,r3)*pow(1-p_err/2,Rd-r1-r2-r3)*p3;

            p_hete = combine_num(Rd,r1)*pow((1-p_err)/2,r1)*pow(1-(1-p_err)/2,Rd-r1)*combine_num(Rd-r1,r2)*pow((1-p_err)/2,r2)*
            pow(1-(1-p_err)/2,Rd-r1-r2)*combine_num(Rd-r1-r2,r3)*pow(p_err,r3)*pow(1-p_err,Rd-r1-r2-r3)*p4;

            if(p_homo>p_hete){
                judge_result[i] = false;
            }
            else{
                judge_result[i] = true;
            }

        }
        else if(dl[i]==4){
            int r1 = cor_seq.base_count[i][0].second;
            int r2 = cor_seq.base_count[i][1].second;
            int r3 = cor_seq.base_count[i][2].second;
            int r4 = cor_seq.base_count[i][3].second;
            p_homo = combine_num(Rd,r1)*pow(1-p_err,r1)*pow(p_err,Rd-r1)*combine_num(Rd-r1,r2)*pow(p_err/3,r2)*pow(1-p_err/3,Rd-r1-r2)*
            combine_num(Rd-r1-r2,r3)*pow(p_err/3,r3)*pow(1-p_err/3,Rd-r1-r2-r3)*combine_num(Rd-r1-r2-r3,r4)*pow(p_err/3,r4)*
            pow(1-p_err/3,Rd-r1-r2-r3-r4)*p3;
            
            p_hete = combine_num(Rd,r1)*pow((1-p_err)/2,r1)*pow(1-(1-p_err)/2,Rd-r1)*combine_num(Rd-r1,r2)*pow((1-p_err)/2,r2)*
            pow(1-(1-p_err)/2,Rd-r1-r2)*combine_num(Rd-r1-r2,r3)*pow(p_err/2,r3)*pow(1-p_err/2,Rd-r1-r2-r3)*combine_num(Rd-r1-r2-r3,r4)*
            pow(p_err/2,r4)*pow(1-p_err/2,Rd-r1-r2-r3-r4)*p4;

            if(p_homo>p_hete){
                judge_result[i] = false;
            }
            else{
                judge_result[i] = true;
            }
        }
        else if(dl[i]==5){
            int r1 = cor_seq.base_count[i][0].second;
            int r2 = cor_seq.base_count[i][1].second;
            int r3 = cor_seq.base_count[i][2].second;
            int r4 = cor_seq.base_count[i][3].second;
            int r5 = cor_seq.base_count[i][4].second;
            p_homo = combine_num(Rd,r1)*pow(1-p_err,r1)*pow(p_err,Rd-r1)*combine_num(Rd-r1,r2)*pow(p_err/4,r2)*pow(1-p_err/4,Rd-r1-r2)*
            combine_num(Rd-r1-r2,r3)*pow(p_err/4,r3)*pow(1-p_err/4,Rd-r1-r2-r3)*combine_num(Rd-r1-r2-r3,r4)*pow(p_err/4,r4)*
            pow(1-p_err/4,Rd-r1-r2-r3-r4)*combine_num(Rd-r1-r2-r3-r4,r5)*pow(p_err/4,r5)*pow(1-p_err/4,Rd-r1-r2-r3-r4-r5)*p3;

            p_hete = combine_num(Rd,r1)*pow((1-p_err)/2,r1)*pow(1-(1-p_err)/2,Rd-r1)*combine_num(Rd-r1,r2)*pow((1-p_err)/2,r2)*
            pow(1-(1-p_err)/2,Rd-r1-r2)*combine_num(Rd-r1-r2,r3)*pow(p_err/3,r3)*pow(1-p_err/3,Rd-r1-r2-r3)*combine_num(Rd-r1-r2-r3,r4)*
            pow(p_err/3,r4)*pow(1-p_err/3,Rd-r1-r2-r3-r4)*combine_num(Rd-r1-r2-r3-r4,r5)*pow(p_err/3,r5)*pow(1-p_err/3,Rd-r1-r2-r3-r4-r5)*p4;

            if(p_homo>p_hete){
                judge_result[i] = false;
            }
            else{
                judge_result[i] = true;
            }
        }
    }
    return judge_result;
}