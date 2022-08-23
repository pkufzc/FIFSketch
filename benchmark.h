#pragma once

#include <chrono>
#include <climits>
#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
#include "util.h"
#include "Elastic.h"
#include "FIFSketch.h"
#include "FIFSketchOpt.h"
#include "CMSketch.h"
#include "CMHeap.h"
#include "CUSketch.h"
#include "CUHeap.h"
#include "HeavyGuardian.h"
#include "USS.h"
#include "SS.h"
#include "ASketch.h"
#include "Salsa.h"
#include "LogLogFilter.h"
#include "FCMSketch.h"
#include "FCMSketchTopK.h"

typedef std::chrono::high_resolution_clock::time_point TP;

inline TP now() { return std::chrono::high_resolution_clock::now(); }

data_type *read_data(const char *PATH, const count_type length,
                     count_type *cnt) {
	data_type *items = new data_type[length];
	data_type *it = items;

	TIMESTAMP *timestamps = new TIMESTAMP[length];
	TIMESTAMP *timestamp = timestamps;

	FILE *data = fopen(PATH, "rb");
	*cnt = 0;
	while (fread(it++, sizeof(data_type), 1, data) > 0 && fread(timestamp++, sizeof(TIMESTAMP), 1, data) > 0) {
		(*cnt)++;
	}

	fclose(data);

	return items;
}

void BenchTopKFlowSize(const char *PATH, uint32_t K) {
	std::cout << "TopK Flow Size Estimation Benchmark" << std::endl << std::endl;
	count_type cnt;
	data_type *items = read_data(PATH, 100000000, &cnt);
	std::cout << "The number of packet:" << cnt << std::endl;
	constexpr int32_t mem_base = 0;
	constexpr int32_t mem_inc = 200000;
	constexpr int32_t mem_var = 5;
	constexpr int32_t cmp_num = 11;
	constexpr int32_t COUNTER_PER_BUCKET = 8;

	Abstract *sketches[mem_var][cmp_num];

	for (int i = 0; i < mem_var; ++i) {
		sketches[i][0] = new Elastic<COUNTER_PER_BUCKET>((i + 1) * mem_inc);
		sketches[i][1] = new FIFSketch<COUNTER_PER_BUCKET>((i + 1) * mem_inc);
		sketches[i][2] = new HeavyGuardian((i + 1) * mem_inc);
		sketches[i][3] = new USS((i + 1) * mem_inc / 100); 
		sketches[i][4] = new SS((i + 1) * mem_inc / 100); 
		sketches[i][5] = new CMHeap((i + 1) * mem_inc);
		sketches[i][6] = new CUHeap((i + 1) * mem_inc);
		sketches[i][7] = new ASketch((i + 1) * mem_inc);
		sketches[i][8] = new SalsaCM((i + 1) * mem_inc);
		sketches[i][9] = new LogLogFilter((i + 1) * mem_inc);
		sketches[i][10] = new FCMSketchTopK<COUNTER_PER_BUCKET>((i + 1) * mem_inc);
	
	}

	// Ground truth
	HashMap mp;
	for (int l = 0; l < cnt; ++l) {
		if (mp.find(items[l]) == mp.end())
			mp[items[l]] = 1;
		else
			mp[items[l]] += 1;
	}
	std::cout << "The number of flow: " << mp.size() << std::endl;
	for (int i = 0; i < mem_var; ++i) {
		int memory = (mem_base + mem_inc * (i + 1)) / 1000;
		std::cout << "Memory size: " << memory
		          << "KB" << std::endl
		          << std::endl;
		for (int l = 0; l < cnt; ++l) {
			for (int j = 0; j < cmp_num; ++j) {
				sketches[i][j]->Insert(items[l]);
			}
		}
		
		ofstream out;
		out.open("./result/TopK_flow_size_estimation_"+to_string(memory)+"KB.txt",ios::out | ios::trunc);
		for (int j = 0; j < cmp_num; ++j) {
			sketches[i][j]->CompareFlowSize(mp, K, out);
			delete sketches[i][j];
		}
		out.close();
		std::cout << std::endl;
	}

	delete items;
}

void BenchAllFlowSize(const char *PATH) {
	std::cout << "All Flow Size Estimation Benchmark" << std::endl << std::endl;
	count_type cnt;
	data_type *items = read_data(PATH, 100000000, &cnt);
	std::cout << "The number of packet:" << cnt << std::endl;
	constexpr int32_t mem_base = 0;
	constexpr int32_t mem_inc = 200000;
	constexpr int32_t mem_var = 5;
	constexpr int32_t cmp_num = 11;
	constexpr int32_t COUNTER_PER_BUCKET = 8;

	Abstract *sketches[mem_var][cmp_num];

	for (int i = 0; i < mem_var; ++i) {
		sketches[i][0] = new Elastic<COUNTER_PER_BUCKET>((i + 1) * mem_inc);
		sketches[i][1] = new FIFSketchOpt<COUNTER_PER_BUCKET>((i + 1) * mem_inc);
		sketches[i][2] = new HeavyGuardian((i + 1) * mem_inc);
		sketches[i][3] = new USS((i + 1) * mem_inc / 100); 
		sketches[i][4] = new SS((i + 1) * mem_inc / 100); 
		sketches[i][5] = new CMSketch((i + 1) * mem_inc);
		sketches[i][6] = new CUSketch((i + 1) * mem_inc);
		sketches[i][7] = new ASketch((i + 1) * mem_inc);
		sketches[i][8] = new SalsaCM((i + 1) * mem_inc);
		sketches[i][9] = new LogLogFilter((i + 1) * mem_inc);
		sketches[i][10] = new FCMSketch((i + 1) * mem_inc);
	
	}

	// Ground truth
	HashMap mp;
	for (int l = 0; l < cnt; ++l) {
		if (mp.find(items[l]) == mp.end())
			mp[items[l]] = 1;
		else
			mp[items[l]] += 1;
	}
	std::cout << "The number of flow: " << mp.size() << std::endl;
	for (int i = 0; i < mem_var; ++i) {
		int memory = (mem_base + mem_inc * (i + 1)) / 1000;
		std::cout << "Memory size: " << memory
		          << "KB" << std::endl
		          << std::endl;
		for (int l = 0; l < cnt; ++l) {
			for (int j = 0; j < cmp_num; ++j) {
				sketches[i][j]->Insert(items[l]);
			}
		}
		
		ofstream out;
		out.open("./result/All_flow_size_estimation_"+to_string(memory)+"KB.txt",ios::out | ios::trunc);
		for (int j = 0; j < cmp_num; ++j) {
			sketches[i][j]->CompareFlowSize(mp, mp.size(),out); 
			delete sketches[i][j];
		}
		out.close();
		std::cout << std::endl;
	}

	delete items;
}

void BenchHH(const char *PATH) {
	std::cout << "Heavy Hitter Benchmark" << std::endl << std::endl;
	count_type cnt;
	data_type *items = read_data(PATH, 100000000, &cnt);
	std::cout << "The number of packet:" << cnt << std::endl;
	constexpr int32_t mem_base = 0;
	constexpr int32_t mem_inc = 200000;
	constexpr int32_t mem_var = 5;
	constexpr int32_t cmp_num = 11;
	constexpr int32_t COUNTER_PER_BUCKET = 8;

	Abstract *sketches[mem_var][cmp_num];

	for (int i = 0; i < mem_var; ++i) {
		sketches[i][0] = new Elastic<COUNTER_PER_BUCKET>((i + 1) * mem_inc);
		sketches[i][1] = new FIFSketch<COUNTER_PER_BUCKET>((i + 1) * mem_inc);
		sketches[i][2] = new HeavyGuardian((i + 1) * mem_inc);
		sketches[i][3] = new USS((i + 1) * mem_inc / 100); 
		sketches[i][4] = new SS((i + 1) * mem_inc / 100); 
		sketches[i][5] = new CMHeap((i + 1) * mem_inc);
		sketches[i][6] = new CUHeap((i + 1) * mem_inc);
		sketches[i][7] = new ASketch((i + 1) * mem_inc);
		sketches[i][8] = new SalsaCM((i + 1) * mem_inc);
		sketches[i][9] = new LogLogFilter((i + 1) * mem_inc);
		sketches[i][10] = new FCMSketchTopK<COUNTER_PER_BUCKET>((i + 1) * mem_inc);
	
	}

	// Ground truth
	HashMap mp;
	for (int l = 0; l < cnt; ++l) {
		if (mp.find(items[l]) == mp.end())
			mp[items[l]] = 1;
		else
			mp[items[l]] += 1;
	}
	std::cout << "The number of flow: " << mp.size() << std::endl;
	for (int i = 0; i < mem_var; ++i) {
		int memory = (mem_base + mem_inc * (i + 1)) / 1000;
		std::cout << "Memory size: " << memory
		          << "KB" << std::endl
		          << std::endl;
		for (int l = 0; l < cnt; ++l) {
			for (int j = 0; j < cmp_num; ++j) {
				sketches[i][j]->Insert(items[l]);
			}
		}
		
		ofstream out;
		out.open("./result/Heavy_hitter_"+to_string(memory)+"KB.txt",ios::out | ios::trunc);
		for (int j = 0; j < cmp_num; ++j) {
			sketches[i][j]->CompareHH(mp, cnt, out);
			delete sketches[i][j];
		}
		out.close();
		std::cout << std::endl;
	}

	delete items;
}

HashMap DiffHC(HashMap map1,HashMap map2){
	HashMap ret;
	for (auto it = map1.begin();it != map1.end();++it){
		ret[it->first] = abs((int)(it->second - map2[it->first]));
	}

	for (auto it = map2.begin();it != map2.end();++it){
		if(ret.find(it->first) == ret.end()){
			ret[it->first] = abs((int)(it->second - map1[it->first]));
		}
	}

	return ret;
}

void CompareHC(HashMap mp, HashMap record, count_type length, ofstream &out) {
	double realHH = 0, estHH = 0, bothHH = 0, aae = 0, are = 0, precision = 0, recall = 0, f1 = 0, alpha = 0.0005;
	count_type threshold = (count_type)(length * alpha);

	for(auto it = mp.begin();it != mp.end();++it){
		bool real, est;
		double realF = it->second, estF = record[it->first];

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

	printf("AAE: %lf\tARE: %lf\tPrecision: %lf\tRecall: %lf\tF1: %lf\n", aae, are, precision, recall, f1);
	out << aae << "\t" << are << "\t" << precision << "\t" << recall << "\t" << f1 << "\n";
}

void BenchHC(const char *PATH) {
	std::cout << "Heavy Change Benchmark" << std::endl << std::endl;
	count_type cnt;
	data_type *items = read_data(PATH, 100000000, &cnt);
	std::cout << "The number of packet:" << cnt << std::endl;
	constexpr int32_t mem_base = 0;
	constexpr int32_t mem_inc = 200000;
	constexpr int32_t mem_var = 5;
	constexpr int32_t cmp_num = 11;
	constexpr int32_t COUNTER_PER_BUCKET = 8;

	Abstract *sketches1[mem_var][cmp_num];
	Abstract *sketches2[mem_var][cmp_num];

	for (int i = 0; i < mem_var; ++i) {
		sketches1[i][0] = new Elastic<COUNTER_PER_BUCKET>((i + 1) * mem_inc);
		sketches1[i][1] = new FIFSketch<COUNTER_PER_BUCKET>((i + 1) * mem_inc);
		sketches1[i][2] = new HeavyGuardian((i + 1) * mem_inc);
		sketches1[i][3] = new USS((i + 1) * mem_inc / 100); 
		sketches1[i][4] = new SS((i + 1) * mem_inc / 100); 
		sketches1[i][5] = new CMHeap((i + 1) * mem_inc);
		sketches1[i][6] = new CUHeap((i + 1) * mem_inc);
		sketches1[i][7] = new ASketch((i + 1) * mem_inc);
		sketches1[i][8] = new SalsaCM((i + 1) * mem_inc);
		sketches1[i][9] = new LogLogFilter((i + 1) * mem_inc);
		sketches1[i][10] = new FCMSketchTopK<COUNTER_PER_BUCKET>((i + 1) * mem_inc);

		sketches2[i][0] = new Elastic<COUNTER_PER_BUCKET>((i + 1) * mem_inc);
		sketches2[i][1] = new FIFSketch<COUNTER_PER_BUCKET>((i + 1) * mem_inc);
		sketches2[i][2] = new HeavyGuardian((i + 1) * mem_inc);
		sketches2[i][3] = new USS((i + 1) * mem_inc / 100); 
		sketches2[i][4] = new SS((i + 1) * mem_inc / 100); 
		sketches2[i][5] = new CMHeap((i + 1) * mem_inc);
		sketches2[i][6] = new CUHeap((i + 1) * mem_inc);
		sketches2[i][7] = new ASketch((i + 1) * mem_inc);
		sketches2[i][8] = new SalsaCM((i + 1) * mem_inc);
		sketches2[i][9] = new LogLogFilter((i + 1) * mem_inc);
		sketches2[i][10] = new FCMSketchTopK<COUNTER_PER_BUCKET>((i + 1) * mem_inc);
	
	}

	HashMap mp1[mem_var][cmp_num],mp2[mem_var][cmp_num],record1[mem_var][cmp_num],record2[mem_var][cmp_num];
	for (int l = 0; l < (cnt >> 1); ++l){
		for (int i = 0; i < mem_var; ++i) {
			for (int j = 0; j < cmp_num; ++j) {
				mp1[i][j][items[l]] += 1;
				sketches1[i][j]->Insert(items[l]);
			}
		}
	}

	for (int l = (cnt >> 1); l < cnt; ++l){
		for (int i = 0; i < mem_var; ++i) {
			for (int j = 0; j < cmp_num; ++j) {
				mp2[i][j][items[l]] += 1;
				sketches2[i][j]->Insert(items[l]);
			}
		}
	}

	for (int i = 0; i < mem_var; ++i) {
		for (int j = 0; j < cmp_num; ++j) {
			for (auto it = mp1[i][j].begin();it != mp1[i][j].end();++it){
				record1[i][j][it->first] = sketches1[i][j]->Query(it->first);
			}
		}
	}

	for (int i = 0; i < mem_var; ++i) {
		for (int j = 0; j < cmp_num; ++j) {
			for (auto it = mp2[i][j].begin();it != mp2[i][j].end();++it){
				record2[i][j][it->first] = sketches2[i][j]->Query(it->first);
			}
		}
	}

	for (int i = 0; i < mem_var; ++i) {
		int memory = (mem_base + mem_inc * (i + 1)) / 1000;
		std::cout << "Memory size: " << memory
		          << "KB" << std::endl
		          << std::endl;
		ofstream out;
		out.open("./result/Heavy_change_"+to_string(memory)+"KB.txt",ios::out | ios::trunc);
		for (int j = 0; j < cmp_num; ++j) {
			std::cout << (sketches1[i][j]->name) << std::endl;
			HashMap realMap = DiffHC(mp1[i][j],mp2[i][j]);
			HashMap estMap = DiffHC(record1[i][j],record2[i][j]);

			CompareHC(realMap, estMap, realMap.size(), out);
		}
		out.close();
	}

	for (int i = 0; i < mem_var; ++i) {
		for (int j = 0; j < cmp_num; ++j) {
			delete sketches1[i][j];
			delete sketches2[i][j];
		}
	}

	delete items;
}

void BenchThp(const char *PATH) {
	std::cout << "Throughput Benchmark"
	          << std::endl
	          << std::endl;

	count_type cnt;
	data_type *items = read_data(PATH, 100000000, &cnt);

	constexpr int32_t mem_base = 0;
	constexpr int32_t mem_inc = 200000;
	constexpr int32_t mem_var = 5;
	constexpr int32_t cmp_num = 15;
	constexpr int32_t round = 5;
	constexpr int32_t COUNTER_PER_BUCKET = 8;

	Abstract *sketches[mem_var][cmp_num];

	int progress = 0;

	for (int i = 0; i < mem_var; ++i) {
		int memory = (mem_base + mem_inc * (i + 1)) / 1000;
		std::cout << "Memory size: " << memory
		          << "KB" << std::endl
		          << std::endl;

		double thp[round][cmp_num] = {};
		double avg_thp[cmp_num] = {};

		for (int j = 0; j < round; ++j) {
			sketches[i][0] = new Elastic<COUNTER_PER_BUCKET>((i + 1) * mem_inc);
			sketches[i][1] = new FIFSketch<COUNTER_PER_BUCKET>((i + 1) * mem_inc);
			sketches[i][2] = new FIFSketchOpt<COUNTER_PER_BUCKET>((i + 1) * mem_inc);
			sketches[i][3] = new HeavyGuardian((i + 1) * mem_inc);
			sketches[i][4] = new USS((i + 1) * mem_inc / 100); 
			sketches[i][5] = new SS((i + 1) * mem_inc / 100); 
			sketches[i][6] = new CMSketch((i + 1) * mem_inc);
			sketches[i][7] = new CUSketch((i + 1) * mem_inc);
			sketches[i][8] = new CMHeap((i + 1) * mem_inc);
			sketches[i][9] = new CUHeap((i + 1) * mem_inc);
			sketches[i][10] = new ASketch((i + 1) * mem_inc);
			sketches[i][11] = new SalsaCM((i + 1) * mem_inc);
			sketches[i][12] = new LogLogFilter((i + 1) * mem_inc);
			sketches[i][13] = new FCMSketch((i + 1) * mem_inc);
			sketches[i][14] = new FCMSketchTopK<COUNTER_PER_BUCKET>((i + 1) * mem_inc);
			

			for (int l = 0; l < cmp_num; ++l) {
				TP start, finish;

				start = now();
				for (int m = 0; m < cnt; ++m) {
					sketches[i][l]->Insert(items[m]);
				}
				finish = now();

				thp[j][l] =
				    (double)cnt /
				    std::chrono::duration_cast< std::chrono::duration<
				        double, std::ratio< 1, 1000000 > > >(finish - start)
				        .count();
				avg_thp[l] += thp[j][l];

				if (j != round - 1) {
					delete sketches[i][l];
				}
			}
		}

		for (int l = 0; l < cmp_num; ++l) {
			printf("%12s:\tthp = %lf\n", sketches[i][l]->name,
			       avg_thp[l] / round);
			delete sketches[i][l];
		}
		std::cout << std::endl;
	}

	delete items;
}
