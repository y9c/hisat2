#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cstring>
#include <map>
#include <vector>

using namespace std;

char fromBase = '\0';
char toBase = '\0';
bool jump = false;

vector<string> splitString(string &inputString, string delimiter)
{
    vector<string> outputVector;
    size_t startPosition = 0;
    size_t endPosition = 0;
    string info;

    while ((endPosition = inputString.find(delimiter, startPosition)) != string::npos) {
        info = inputString.substr(startPosition, endPosition-startPosition);
        outputVector.push_back(info);
        startPosition = endPosition + delimiter.length();
    }
    outputVector.push_back(inputString.substr(startPosition));
    return outputVector;
}

string joinVector(vector<string> &inputVector, string delimiter)
{
    string outputString = "";
    for (int i = 0; i < inputVector.size()-1; i++){
        outputString += inputVector[i] + delimiter;
    }


    outputString += inputVector[inputVector.size()-1];

    return outputString;
}

char complement[] = {
        /*   0 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /*  16 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /*  32 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,'-',  0,  0,
        /*  48 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /*  64 */ 0,'T','V','G','H',  0,  0,'C','D',  0,  0,'M',  0,'K','N',  0,
        /*           A   B   C   D           G   H           K       M   N */
        /*  80 */ 0,  0,'Y','S','A',  0,'B','W',  0,'R',  0,  0,  0,  0,  0,  0,
        /*               R   S   T       V   W       Y */
        /*  96 */ 0,'T','V','G','H',  0,  0,'C','D',  0,  0,'M',  0,'K','N',  0,
        /*           a   b   c   d           g   h           k       m   n */
        /* 112 */ 0,  0,'Y','S','A',  0,'B','W',  0,'R',  0,  0,  0,  0,  0,  0,
        /*               r   s   t       v   w       y */
        /* 128 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 144 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 160 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 176 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 192 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 208 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 224 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 240 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
};

char convertPlanA[] = {
        /*   0 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /*  16 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /*  32 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,'-',  0,  0,
        /*  48 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /*  64 */ 0,'A','B','C','D',  0,  0,'G','H',  0,  0,'K',  0,'M','N',  0,
        /*           A   B   C   D           G   H           K       M   N */
        /*  80 */ 0,  0,'R','S','T',  0,'V','W',  0,'Y',  0,  0,  0,  0,  0,  0,
        /*               R   S   T       V   W       Y */
        /*  96 */ 0,'A','B','C','D',  0,  0,'G','H',  0,  0,'K',  0,'M','N',  0,
        /*           a   b   c   d           g   h           k       m   n */
        /* 112 */ 0,  0,'R','S','T',  0,'V','W',  0,'Y',  0,  0,  0,  0,  0,  0,
        /*               r   s   t       v   w       y */
        /* 128 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 144 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 160 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 176 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 192 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 208 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 224 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 240 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
};

char convertPlanB[] = {
        /*   0 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /*  16 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /*  32 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,'-',  0,  0,
        /*  48 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /*  64 */ 0,'A','B','C','D',  0,  0,'G','H',  0,  0,'K',  0,'M','N',  0,
        /*           A   B   C   D           G   H           K       M   N */
        /*  80 */ 0,  0,'R','S','T',  0,'V','W',  0,'Y',  0,  0,  0,  0,  0,  0,
        /*               R   S   T       V   W       Y */
        /*  96 */ 0,'A','B','C','D',  0,  0,'G','H',  0,  0,'K',  0,'M','N',  0,
        /*           a   b   c   d           g   h           k       m   n */
        /* 112 */ 0,  0,'R','S','T',  0,'V','W',  0,'Y',  0,  0,  0,  0,  0,  0,
        /*               r   s   t       v   w       y */
        /* 128 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 144 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 160 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 176 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 192 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 208 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 224 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 240 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
};

void changeBase (string &str, char changedStrA[], char changedStrB[])
{
    //changedStrA.clear();
    //changedStrB.clear();

    for (int i = 0; i < str.size(); i++) {
        changedStrA[i] = convertPlanA[str[i]];
        changedStrB[i] = convertPlanB[str[i]];
    }

    changedStrA[str.size()] = '\0';
    changedStrB[str.size()] = '\0';
}


void argumentParser (int argc, const char *argv[], string & inFileName, string & outFileName1, string & outFileName2){

    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "--reference") == 0){
            inFileName = argv[i+1];
            continue;
        }
        if (strcmp(argv[i], "--o1") == 0){
            outFileName1 = argv[i+1];
            continue;
        }
        if (strcmp(argv[i], "--o2") == 0){
            outFileName2 = argv[i+1];
            continue;
        }
        if (strcmp(argv[i], "--base-change") == 0){
            string baseChange = argv[i+1];
            fromBase = baseChange.front();
            toBase = baseChange.back();
            continue;
        }
    }

    convertPlanA[fromBase] = toBase;
    convertPlanA[tolower(fromBase)] = toBase;
    convertPlanB[complement[fromBase]] = complement[toBase];
    convertPlanB[tolower(complement[fromBase])] = complement[toBase];
}


void changeBaseInFasta(string &inFileName, ofstream &outFile1, ofstream &outFile2)
{

    ifstream inFile;
    inFile.open(inFileName, ios_base::in);

    if (inFile.is_open())
    {
        string line;
        vector<string> lineVector;

        char changedLine1[1000];
        char changedLine2[1000];

        while (inFile.good())
        {
            getline(inFile, line);
            if (line.front() != '>'){
                changeBase(line, changedLine1, changedLine2);

                if (outFile1.is_open()) {
                    if (!jump){
                        outFile1 << changedLine1 << "\n";
                    }
                }
                if (outFile2.is_open()) {
                    if (!jump){
                        outFile2 << changedLine2 << "\n";
                    }
                }
            }
            else{
                lineVector = splitString(line, " ");
                if (lineVector[0] != ">na"){
                    //lineVector[0] = lineVector[0] + "_" + fromBase + toBase;
                    string title1 = joinVector(lineVector, " ");
                    string title2 = joinVector(lineVector, " ");
                    if (outFile1.is_open()) {
                            outFile1 << title1 << "\n";
                    }
                    if (outFile2.is_open()) {
                            outFile2 << title2 << "\n";
                    }
                    jump = false;
                } else {
                    jump = true;
                }

            }


        }
    }

    inFile.close();
    //outFile.close();
}

int main(int argc, const char * argv[]) {

    string inFileName;
    string outFileName1;
    string outFileName2;

    argumentParser(argc, argv, inFileName, outFileName1, outFileName2);

    ofstream outFile1;
    outFile1.open(outFileName1, ios_base::out);

    ofstream outFile2;
    outFile2.open(outFileName2, ios_base::out);

    //cout << convertPlanA['T'];
    changeBaseInFasta(inFileName, outFile1, outFile2);

    //fromBase = complement[fromBase];
    //toBase = complement[toBase];

    //changeBaseInFasta(inFileName, outFile1, outFile2);

    outFile1.close();
    outFile2.close();

    return 0;
}
