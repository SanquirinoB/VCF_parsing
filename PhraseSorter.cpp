#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>
// #include <stxxl/io>
// #include <stxxl/vector>
// #include <stxxl/stream>

using namespace std;

struct phrase
{
    typedef unsigned char three_digits; // 1 byte
    typedef unsigned short four_digits; // 2 byte
    typedef unsigned ten_digits; // 4 bytes

    // Optimal aligment
    four_digits indv;
    three_digits chrom, alele;
    ten_digits pos;
    ten_digits pos_e;
    char *edit[4], *len[6], *len_e[6];
    bool is_short;

    phrase() { }
    phrase(string line)
    {
        char* buffer[10];
        string aux;
        int int_aux;
        // indv
        //memcpy(buffer, line, 4);
        //aux = &buffer;
        //indv = stoi(aux, nullptr);
        // chrom
        //memcpy(three_d, line + 4, 3);
        // alele
        //memcpy(three_d, line + 4, 3);
        // pos
        // len
        // edit
        // pose_e
        // len_e
    }
};

int main()
{
    cout << alignof(phrase) << endl; // output: 16
}