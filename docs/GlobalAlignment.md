##Approach

Pretty straight forward.
Assuming M[i][j] represents the optimal score for 
aligning two sequences. Then M[i][j]
is given by the famous Needleman Wunsch Dynamic programming
algorithm:

M[i][j] = max(M_<sub>mm</sub>, M___<sub>_i1</sub>>, M_<sub>i2</sub>>)

where M_mm = M[i-1][j-1] + W(U_<sub>i</sub>, U___<sub>_j</sub>>) 

M_  <sub>i1 = M[i-1][j] + W(U_i, -) => Aligning U(U__<sub>i</sub>) against a gap 

M_i1 = M[i-1][j] + W(U-, V__<sub>j</sub> )  =>  Aligning V(V__<sub>j</sub>) against a gap 


## Initialization

We add “_ _ ” before each sequence
