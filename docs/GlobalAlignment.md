##Approach

Pretty straight forward.
Assuming M<sub>ij</sub> represents the optimal score for 
aligning two sequences U:u<sub>1</sub>u<sub>2</sub>....u<sub>m</sub> with 
sequence V: v<sub>1</sub>v<sub>2</sub>....v<sub>n</sub>.
Let M<sub>ij</sub> represent the optimal score of aligning the i<sup>th</sup> position of U
against the j<sup>th</sup> position of V.
Then M<sub>ij</sub> is given by the famous Needleman Wunsch Dynamic programming
algorithm:

M<sub>ij</sub> = max(M<sub>mm</sub>, M<sub>i1</sub>, M<sub>i2</sub>)

where M<sub> = M<sub>(i-1)(j-1)</sub? + W(U<sub>i</sub>, U<sub>j</sub>) 

M<sub>del1</sub> = M<sub>i(j-1)</sub> + W(U<sub>i</sub>, -) => Aligning U(U<sub>i</sub>) against a gap 

M<sub>del2</sub> = M<sub>(i-1)j</sub> + W(-, V<sub>j</sub>)  =>  Aligning V(V<sub>j</sub>) against a gap 

Where W(X,Y) is itself a scoring function given by:

W(X,Y) = 2 if X==Y

W(X,Y) = -1 if X!=Y

W(X,Y) = -2 if X== '-' or Y == '-'


