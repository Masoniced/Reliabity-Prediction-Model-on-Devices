// ConsoleApplication1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
using namespace std;

bool check(string const& s) {
    string first = s.substr(0, 1);
    if (first == "0") {
        if (s.size() == 1) {
            return true;
        }
        else {
            return false;
        }
    }
    else {
        int value = stoi(s);
        if (value < 0 || value > 255) {
            return false;
        }
        else {
            return true;
        }
    }
}

vector<string> output(vector<string> s) {
    vector<string> result;
    int first = stoi(s[0]);
    int len1, len2;
    if (first > 99) {
        len1 = 3;
    }
    else if (first > 9) {
        len1 = 2;
    }
    else {
        len1 = 1;
    }

    int last = stoi(s[3]);
    if (last > 99) {
        len2 = 3;
    }
    else if (last > 9) {
        len2 = 2;
    }
    else {
        len2= 1;
    }

    for (int i = 0; i < len1; i++) {
        for (int j = 0; j < len2; j++) {
            string word1 = s[0].substr(i, 3 - i);
            string word2 = s[3].substr(0, j + 1);
            word1.append(".");
            word1.append(s[1]);
            word1.append(".");
            word1.append(s[2]);
            word1.append(".");
            word1.append(word2);
            result.push_back(word1);
        }
    }
    return result;
}

int main()
{
    string input = "1.1.1111.16..172.16.254.1.1";
    int length = input.size();
    vector<string> fragment;
    int startIndex = 0;
    int endIndex = 0;
    string iterate;
    for (int i = 0; i < length; i++) {
        string temp = input.substr(i, 1);
        if (temp == ".") {
            if (startIndex != endIndex) {
                iterate = input.substr(startIndex, endIndex - startIndex);
            }
            else {
                string empty;
                iterate = empty;
            }
            fragment.push_back(iterate);
            if (i < length - 1) {
                startIndex = i + 1;
                endIndex = i + 1;
            }
        }
        else {
            if (i < length - 1) {
                endIndex += 1;
            }
            else {
                iterate = input.substr(startIndex, endIndex);
                fragment.push_back(iterate);
            }
        }
    }

    int sectionLength = fragment.size();
    vector<int> checkList(sectionLength, 0);
    for (int i = 0; i < sectionLength; i++) {
        if (fragment[i].size() != 0) {
            if (check(fragment[i])) {
                checkList[i] = 1;
            }
        }
    }
    vector<string> result;
    for (int i = 0; i < checkList.size() - 3; i++) {

        if (checkList[i + 1] == 1 && checkList[i + 2] == 1) {
            if (checkList[i] == 1 && checkList[i + 3] == 1) {
                vector<string> word;
                for (int j = 0; j < 4; j++) {
                    word.push_back(fragment[i + j]);
                }
                vector<string> temp = output(word);
                result.insert(result.end(), temp.begin(), temp.end());
            }

            if (checkList[i] == 0) {
                if (fragment[i].size() > 0) {
                    string word1;
                    int first = stoi(fragment[i]);
                    if (first == 0) {
                        word1 = to_string(0);
                    }
                    else {
                        if (first > 255) {
                            word1 = to_string(first % 100);
                        }
                        else {
                            word1 = to_string(first);
                        }
                    }
                    if (checkList[i + 3] == 1) {
                        vector<string> word;
                        word.push_back(word1);
                        for (int j = 1; j < 4; j++) {
                            word.push_back(fragment[i + j]);
                        }
                        vector<string> temp = output(word);
                        result.insert(result.end(), temp.begin(), temp.end());
                    }
                    else {
                        string word2;
                        vector<string> word;
                        if (fragment[i + 3].substr(0, 1) == "0") {
                            word2 = to_string(0);
                        }
                        else {
                            int last = stoi(fragment[i + 3]);
                            if (last < 255) {
                                word2 = to_string(last);
                            }
                            else if (last >= 255 && last < 1000) {
                                word2 = fragment[i+3].substr(0, 2);
                            }
                            else {
                                word2 = fragment[i + 3].substr(0, 3);
                            }
                        }
                        word.push_back(word1);
                        word.push_back(fragment[i + 1]);
                        word.push_back(fragment[i + 2]);
                        word.push_back(word2);
                        vector<string> temp = output(word);
                        result.insert(result.end(), temp.begin(), temp.end());
                    }
                }
            }
        }
    }

    for (int i = 0; i < result.size(); i++) {
        cout << result[i] << endl;
    }

    return 0;
}
