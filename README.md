# OneSketch: A Generic and Accuracy Sketch for Data Streams
## Introduction
In this paper, we aim to propose a generic sketch algorithm capable of achieving more accurately the following five important tasks: finding top-$k$ frequent items, finding heavy hitters, per-item frequency estimation, and heavy changes in time and spatial dimension. State-of-the-art (SOTA) sketch solution for multiple measurement tasks is ElasticSketch (ES), but the accuracy of its frequency estimation still has room for improvement. The reason is that ES introduces errors when querying frequent items, and suffers from overestimation errors for infrequent items. To address these problems, we propose a generic OneSketch for five tasks, with the key design philosophy of one-sided approaching to minimise overestimation errors. Closely surrounding this methodology, we propose four key techniques, which embrace hash collisions and minimize possible errors by well handling extremely recurrent item replacements. Experimental results show that OneSketch clearly outperforms 10 SOTA schemes. For example, compared with ES, OneSketch achieves more than 10$\times$ lower Average Absolute Error on finding top-$k$ frequent items and heavy hitters, and 48.3\% higher F1 Score on heavy changes.  

## Requirements
* cmake
* g++

## How to run
1. cmake .
2. make
3. ./bench
