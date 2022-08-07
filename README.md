# FIFSketch: A Frequent-Item-Friendly Generic Sketch
## Introduction
Finding top-$k$ frequent items in data streams has been a well-studied topic, and sketches are the most acceptable solution for this task.
%
The challenge of handling this task in real time is the trade-off between accuracy and speed. 
%
State-of-the-art sketches can be classified into two categories: record-all-evict-infrequent and separate-and-guard-frequent, and ElasticSketch (ES) is an excellent representative of the latter. 
%
Our proposed FIFSketch aims to enhance ES: it achieves more accurate top-$k$ items frequency estimation, while also more accurate for per-item frequency estimation, heavy hitter and heavy change detection tasks.
%
Experimental results show that FIFSketch obviously outperforms 9 state-of-the-art schemes, $e.g.$, FIFSketch achieves a 13.2\% higher F1 Score than ES under 200KB of memory in top-$k$ items frequency estimation, and an average 5.2 times lower ARE than ES in per-item frequency estimation.
