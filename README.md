# TxRx
Example code for Tx/Rx paper (INSERT TITLE, JOURNAL, AND DOI). 

Two scripts are included:

1) root_placement_strategy.m
  An M=10 element array is studied. The desired two-way response is a -30 dB Chebyshev window. This is achieved via the alternating root placement strategy outlined in the paper. A brute-force method that tests all possible root combinations confirms this result.
2) recursive_formula.m
  An implementation of Eq.(10) in the paper, which is a recursive formula for obtaining a two-way uniform window.
