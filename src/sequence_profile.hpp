#include <ctime>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>

using std::string;
using std::vector;
const char DNA[] = {'A', 'C', 'G', 'T', '_'};
const char AA[] = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'W', 'Y', 'V', '_'};

class ProfileMatrix{
    //Rows correspnoding to 4+1 DNA bases
    //or 20+1 AA.
    //Columns correspond to individual sequence frequency
    int rows, columns;
    //Pointer to a pointer, since an array <=> pointer
    float **P;
    public:
        void initProfileMatrix(){
            // Pointer to pointer p is assigned to an array of integers
            P = new float*[rows];
            for (int i=0; i<rows; i++){
                P[i] = new float [columns];
                for(int j=0; j<columns; j++){
                    P[i][j] = 0.0;
                }
            }
        }
        ProfileMatrix(int rows, int columns);
        ~ProfileMatrix();

};

