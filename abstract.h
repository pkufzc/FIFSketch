#pragma once

#include "hash.h"
#include "util.h"
#include <algorithm>
#include <string>
#include <random>
#include <unordered_map>
#include <climits>
#include <cstring>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <fstream>

using namespace std;

class Abstract {
public:
	char *name;
	struct node {
		data_type x; 
		count_type y;
	};
	node actual[50000000],predicted[50000000];

	Abstract(char *_name) : name(_name){};
	virtual ~Abstract(){};

	virtual void Insert(const data_type item) = 0;
	virtual count_type Query(const data_type item) = 0;

	void CompareFlowSize(HashMap mp, uint32_t k, ofstream &out) {
		
		HashMap actual_order, predicted_order, actual_topk, predicted_topk, all_order;
		count_type hit = 0, cnt1 = 0, cnt2 = 0;
		double aae = 0, are = 0, arre = 0, f1 = 0;
		
		for (HashMap::iterator it = mp.begin(); it != mp.end(); ++it) {
			actual[++cnt1].x = it->first;
        	actual[cnt1].y = it->second;
			predicted[++cnt2].x = it->first;
			predicted[cnt2].y = Query(it->first);
		}

		std::sort(actual+1,actual+cnt1+1,[](node i,node j) {return i.y>j.y;});
    	std::sort(predicted+1,predicted+cnt2+1,[](node i,node j) {return i.y>j.y;});

		for(uint32_t i = 1; i <= cnt1; i++){
			all_order[actual[i].x] = i;
		}

		for(uint32_t i = 1; i <= k; i++){
			actual_order[actual[i].x] = i;
			predicted_order[predicted[i].x] = i;
			actual_topk[actual[i].x] = actual[i].y;
			predicted_topk[predicted[i].x] = predicted[i].y;
		}
		
		for (HashMap::iterator it = actual_topk.begin(); it != actual_topk.end(); ++it) {
			data_type key = it->first;
			count_type value = it->second; 
			if(predicted_topk[key]){
				hit++;
				aae += abs(int(Query(key) - value));
				are += abs(int(Query(key) - value)) / (double)value;
				arre += abs(int(all_order[key] - predicted_order[key])) / (double)all_order[key];
			} 
			
		}

		arre /= hit;
		aae /= hit;
		are /= hit;
		f1 = hit / (double)k;

		printf("%12s:\tAAE: %lf\tARE: %lf\tF1: %lf\n", name, aae, are, f1);
		out << aae << "\t" << are << "\t" << f1 << "\n";

	}

	void CompareHH(HashMap mp, count_type length, ofstream &out) {
		double realHH = 0, estHH = 0, bothHH = 0, aae = 0, are = 0, precision = 0, recall = 0, f1 = 0, alpha = 0.00002;
		count_type threshold = (count_type)(length * alpha);

        for(auto it = mp.begin();it != mp.end();++it){
            bool real, est;
            double realF = it->second, estF = Query(it->first);

            real = (realF > threshold);
            est = (estF > threshold);

            realHH += real;
            estHH += est;

            if(real && est){
                bothHH += 1;
                aae += abs(realF - estF);
                are += abs(realF - estF) / realF;
            }
        }

		aae /= bothHH;
		are /= bothHH;
		precision = bothHH / estHH;
		recall = bothHH / realHH;
		f1 = (2*precision*recall) / (precision + recall);

		printf("%12s:\tAAE: %lf\tARE: %lf\tPrecision: %lf\tRecall: %lf\tF1: %lf\n", name, aae, are, precision, recall, f1);
		out << aae << "\t" << are << "\t" << precision << "\t" << recall << "\t" << f1 << "\n";
	}


};