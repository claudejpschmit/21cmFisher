#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
using namespace std;

//Structure required for sorting according to the first element in the lines.
//Gives the comparison rule for the std::sort algorithm.
struct FirstColumnCmp
{
    bool operator()(const vector<double>& lhs,\
            const vector<double>& rhs) const
    {
        return lhs[0] < rhs[0];
    }
};

int main(int argc, char* argv[])
{
    string filename_in, filename_out;
    if (argc == 3) {
        filename_in = argv[1];
        filename_out = argv[2];
    } else if (argc == 2) {
        filename_in = argv[1];
        filename_out = filename_in; 
    } else {
        cout << "Not the right number of arguments in SortFiles." << endl; 
    }

    vector<vector<double>> file_contents;
    ifstream inFile(filename_in);
    string line;
    while (getline(inFile, line))
    {
        istringstream iss(line);
        int l;
        double Fl, cond_n;
        // This tests whether the line contains exactly 3 elements.
        // It only reads the data if it does, ie. ignoring empty lines.
        if ((iss >> l >> Fl >> cond_n))
        {
            vector<double> row;
            row.push_back(l);
            row.push_back(Fl);
            row.push_back(cond_n);

            file_contents.push_back(row);
        }
    }
    inFile.close();

    //Sorting by std::sort in <algorithm>
    sort(file_contents.begin(), file_contents.end(), FirstColumnCmp());

    ofstream outFile(filename_out);
    for (int i = 0; i < file_contents.size(); i++)
    {
        outFile << file_contents[i][0] << " " <<\
            file_contents[i][1] << " " <<\
            file_contents[i][2] << endl;
    }
    outFile.close();
    return 0;
}
