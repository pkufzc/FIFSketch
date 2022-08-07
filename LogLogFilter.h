#pragma once

#include <iostream>
#include "abstract.h"
#include "CUSketch.h"
#include "LogFilter.h"

class LogLogFilter : public Abstract {
public:
	LogLogFilter(double _Memory) : Abstract((char *)"LogLogFilter") {
		cu = new CUSketch(_Memory * CU_RATIO);
		lf = new LogFilter(_Memory * LF_RATIO);
	}

	~LogLogFilter() {
		delete cu;
		delete lf;
	}
    
	void Insert(const data_type key)
	{
		if (lf->probe(key)) {
			lf->insert(key);
		}
		else {
			cu->Insert(key);
		}
		lf->timestamp++;
	}

	count_type Query(const data_type key) {
		double ret = 0.0;
		if (lf->probe(key)) {
			ret = lf->query1(key);
		}
		else {
			ret = lf->query2() + (double)cu->Query(key);
		}			
		return (uint32_t)ret;
	}
    
private:
	const double CU_RATIO = 0.75;
	const double LF_RATIO = 0.25;
	CUSketch *cu;
	LogFilter *lf;
};
