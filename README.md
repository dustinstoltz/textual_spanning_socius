# Textual Spanning reproduction guide
===========================

This is code and data necessary to reproduce the graphs for Stoltz and Taylor (2019) "Textual Spanning: Finding Discursive Holes in Text Networks" in _Socius_.

The measure of textual spanning we propose works on a document by document similarity matrix. The basic data structure in text analysis is an MxN matrices of documents by terms, n-grams, parts of speech, topics and so on. The rows, therefore, are vector representations of the text and can be easily compared for similarity, usually with cosine similarity. The result is a one-mode document by document matrix, which can easily be interpreted as a weighted adjacency matrix amendable to network metrics. 

The `textSpan` function takes this document by document similarity matrix and outputs a document specific measure which increases when a document is similar to documents which are not also similar to each other. This is defined by the following equations:
<img src="https://latex.codecogs.com/gif.latex?S_i%20%3D%20%5Csum_j%20%5Cleft%20%28%20p_%7Bij%7D%20&plus;%20%5Csum_q%20%5Cfrac%7Bp_%7Bqj%7D%7D%7Bp_%7Biq%7D%7D%20%5Cright%20%29%5E2"/>
