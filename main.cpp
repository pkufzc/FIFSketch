#include "benchmark.h"
#include <string>
#include <iostream>

using namespace std;

const string folder = "./datasets/";
const string filenames[] = {"130000.dat"};

int main() {
	cout << endl << "**Benchmark**" << endl << endl;
	uint32_t K;
	cin >> K;
	BenchTopKFlowSize((folder + filenames[0]).c_str(), K);
	BenchAllFlowSize((folder + filenames[0]).c_str());
	BenchHH((folder + filenames[0]).c_str());
	BenchHC((folder + filenames[0]).c_str());
	BenchThp((folder + filenames[0]).c_str());
}