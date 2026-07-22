#ifndef SIMPLETOOLS_H
#define SIMPLETOOLS_H

#include <cmath>


using namespace std;

TString print_vec(vector<TString> vec) {
    TString out = "";
    for (int i = 0; i < vec.size(); i++) {
        if (i == vec.size()-1) out += vec[i] + ".";
        else out += vec[i] + ", ";
    }
    return out;
}

double median(vector<double> v) {
    sort(v.begin(), v.end());
    if (v.size()%2 == 1) {
        return v[v.size()/2];
    }
    else {
        return .5 * (v[v.size()/2] + v[(v.size()/2)-1]);
    }
}

bool contains_TString(vector<TString> vec, TString x) {
    return (std::find(vec.begin(), vec.end(), x) != vec.end());
}

// return a TString representation of a double, rounded to nearest .1
TString to_string_round(double num) {
    if (num < 1) {
        return "." + to_string((int)(10*num));
    }
    else {
        return to_string((int)num) + "." + to_string((int)(10*fmod(num, 1.0)));
    }
}

void print_bin (std::vector<int> bin) {
    cout << "pt of " << bin[0] << " to " << bin[1] << ", " << bin[2] << " events." << endl;
}

void print_bins (std::vector<vector<int>> bins) {
    for (int iBin = 0; iBin < bins.size(); iBin++) {
        cout << "Bin " << iBin << ": ";
        print_bin(bins[iBin]);
    }
}

void print_bins2code (std::vector<vector<int>> bins) {
    cout << "vector<vector<int>> bins = {" << endl;
    for (int iBin = 0; iBin < bins.size(); iBin++) {
        cout << "{" << bins[iBin][0] << ", " << bins[iBin][1] << ", " << bins[iBin][2] << (iBin == bins.size()-1? "}" : "},") << endl;
    }
    cout << "};" << endl;

}


TString bool2Str(bool b) {return b? "True" : "False";}



#endif