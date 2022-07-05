/*
 * Use byte code for represenation of sequences. or two boolean
 *
 */

#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

static const int N = 10000;
static const int M = N;
static const bool SHOWMATRIX = false;
static const bool SHOWALIGNMENT = true;
static const int MATCH = 1;
static const int MIS = -1;
static const int GAP = -2;
static const int BENCHMARK = 0;
static const char DNACHARS[] = "ATCG";
static const int LENALPHABET = sizeof(DNACHARS);
static const int MSA[3][2] = {{0, 1}, {0, 2}, {1, 2}};
static const string DATALOC = "DATA/sequences.txt";

class Align
{
    int row;
    int col;

    char* ps1;
    char* ps2;
    std::vector<int> matrix;
    std::vector<char> tmatrix;

public:
    Align(int m, int n, bool verboseAlignment = false, bool verboseMatrix = false, const int *MSA = NULL) : 
        row(m), col(n), ps1(new char[n]), ps2(new char[m]), 
        matrix((n+1)*(m+1)), tmatrix((n+1)*(m+1))
    {
        if (!MSA)
        {
            this->_getRandomDNA(ps1, col);
            this->_getRandomDNA(ps2, row);
        }
        else
        {
            ps1 = this->_readSequence(ps1, MSA[0]);
            ps2 = this->_readSequence(ps2, MSA[1]);
        }
        row++;
        col++;

        this->_fillMatrix(matrix.data(), tmatrix.data(), ps1, ps2);

        if (verboseMatrix)
            this->_verboseMatrix(matrix.data(), tmatrix.data(), ps1, ps2);
        
        this->_traceback(tmatrix.data(), ps1, ps2);

        if (verboseAlignment){
            this->_verboseAlignment();
            cout << "Alignment Score: " << matrix[row * col - 1] << endl;
        }
    };

    

    ~Align(){
        delete [] ps1;
        delete [] ps2;
    }

private:
    static void _getRandomDNA(char strand[], int size)
    {
        for (unsigned int i = 0; i < size; i++)
            strand[i] = DNACHARS[rand() % (LENALPHABET - 1)];
    };

    static char *_readSequence(char strand[], int size)
    {
        int line_no = 0;
        ifstream sFile;
        sFile.open(DATALOC);
        if (!sFile.is_open())
            throw std::runtime_error("Can't open file: " + DATALOC);
        string sLine;
        while (line_no != (size + 4))
        {
            ++line_no;
            getline(sFile, sLine);
        }
        strcpy(strand, sLine.c_str());
        return strand;
    };

    void _verboseMatrix(int mat[], char tmat[], char s1[], char s2[])
    {
        this->_printMatrix(mat, s1, s2);
        cout << endl
             << endl;
        this->_printMatrix(tmat, s1, s2);
        cout << endl;
    };

    void _verboseAlignment()
    {

        ifstream s1File, s2File;
        s1File.open("s1.txt");
        s2File.open("s2.txt");
        string line1, line2;
        int c = 0;
        for (; getline(s1File, line1) && getline(s2File, line2);)
        {
            cout << "query " << (c + 1) << " " << line1
                 << " " << (c + 80) << endl
                 << "sbjct " << (c + 1) << " "
                 << line2 << " " << (c + 80) << endl
                 << endl;
            c += 80;
        }
    };

    int *_fillMatrix(int mat[], char tmat[], char s1[], char s2[])
    {
        int index;
        for (unsigned int i = 0; i < col; i++)
        {
            index = (i) + col * (0);
            mat[index] = GAP * i;
            tmat[index] = 'H';
        }

        for (unsigned int i = 0; i < row; i++)
        {
            index = (0) + col * (i);
            mat[index] = GAP * i;
            tmat[index] = 'V';
        }
        int ma, ii, jj, idx;
        int c[3];

        const static char s[3] = {'D', 'H', 'V'};

        for (unsigned int j = 1; j < row; j++)
        {
            for (unsigned int i = 1; i < col; i++)
            {
                ii = i-1; jj = j-1;
                ma = s1[ii] == s2[jj];

                c[0] /*di*/ = mat[(ii) + col * (jj)] + ((ma * MATCH) + ((ma-1) * (-MIS)));
                c[1] /*ho*/ = mat[(ii) + col * (j)] + GAP;
                c[2] /*ve*/ = mat[i + col * (jj)] + GAP;

                index = i + col * j;
                
                idx = max_element(c, c + 3) - c;
                mat[index] = c[idx];
                tmat[index] = s[idx];
            }
        }

        return mat;
    };

    void _traceback(char tmat[], char s1[], char s2[])
    {
        ofstream s1File, s2File;
        s1File.open("s1.txt");
        s2File.open("s2.txt");
        tmat[0] = 'e';
        int i = 1;
        int j = 1;
        char current;
        int nChars = 0;
        do
        {
            current = tmat[(col - i) + col * (row - j)];
            nChars++;
            if (current == 'D')
            {
                i++;
                j++;
                s1File << s1[col - i];
                s2File << s2[row - j];
            }
            else if (current == 'V')
            {
                j++;
                s1File << "-";
                s2File << s2[row - i];
            }
            else if (current == 'H')
            {
                i++;
                s1File << s1[col - i];
                s2File << "-";
            }
            if (nChars % 80 == 0)
            {
                s1File << endl;
                s2File << endl;
            }

        } while (current != 'e');

        s1File.close();
        s2File.close();
    }
    
    void _printMatrix(auto mat[], char s1[], char s2[])
    {
        cout << "      ";
        for (unsigned int i = 0; i < col - 1; i++)
        {
            cout << s1[i] << ", ";
        }
        char n;
        for (unsigned int i = 0; i < row; i++)
        {
            n = i != 0 ? s2[i - 1] : ' ';
            cout << endl
                 << " " << n << " ";
            for (unsigned int j = 0; j < col; j++)
            {
                cout << mat[(j) + col * (i)] << ",";
            }
        }
    };
};

void benchMark(int times)
{
    double sum = 0.0;
    for (unsigned int i = 0; i < times; i++)
    {
        clock_t begin = clock();

        Align(N, M);

        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        sum += elapsed_secs;
    }

    cout << "Average elapsed time(s) "
         << times << " iterations: "
         << (sum / times) << endl;
}

int main()
{
    srand(1);
    if (BENCHMARK != 0)
    {
        benchMark(BENCHMARK);
    }
    else
    {
        clock_t begin = clock();
        for (int i = 0; i < 1; i++)
        {
            Align(N, M, SHOWALIGNMENT, SHOWMATRIX, NULL);
        }
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        cout << "Elapsed time(s) " << elapsed_secs << endl;
    }
}
