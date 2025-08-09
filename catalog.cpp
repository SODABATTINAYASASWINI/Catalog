#include <bits/stdc++.h>
#include "nlohmann/json.hpp"
using namespace std;
using json = nlohmann::json;

// Decode value string in given base (up to 36) to long long
long long decodeValue(const string& val, int base) {
    long long res = 0;
    for (char ch : val) {
        int digit;
        if (isdigit(ch)) digit = ch - '0';
        else digit = tolower(ch) - 'a' + 10;
        res = res * base + digit;
    }
    return res;
}

// Solve least squares for A*x= -C with c=1 fixed, returns a,b
pair<double,double> solveForAandB(const vector<vector<double>>& A, const vector<double>& C) {
    int n = A.size();
    vector<vector<double>> mat(n, vector<double>(2));
    vector<double> rhs(n);
    for (int i = 0; i < n; i++) {
        mat[i][0] = A[i][0]; // coefficient of a (x^2)
        mat[i][1] = A[i][1]; // coefficient of b (y)
        rhs[i] = -C[i];      // -c, c=1 here
    }

    // Compute ATA and ATb
    vector<vector<double>> ATA(2, vector<double>(2, 0));
    vector<double> ATb(2, 0);
    for (int i = 0; i < n; i++) {
        ATA[0][0] += mat[i][0] * mat[i][0];
        ATA[0][1] += mat[i][0] * mat[i][1];
        ATA[1][0] += mat[i][1] * mat[i][0];
        ATA[1][1] += mat[i][1] * mat[i][1];

        ATb[0] += mat[i][0] * rhs[i];
        ATb[1] += mat[i][1] * rhs[i];
    }

    double det = ATA[0][0] * ATA[1][1] - ATA[0][1] * ATA[1][0];
    if (fabs(det) < 1e-12) return {0, 0}; // singular matrix fallback

    double a = (ATb[0] * ATA[1][1] - ATb[1] * ATA[0][1]) / det;
    double b = (ATA[0][0] * ATb[1] - ATA[1][0] * ATb[0]) / det;

    return {a, b};
}

// Process one test case JSON object
void processTestCase(const json& test) {
    int n = test["keys"]["n"];
    int k = test["keys"]["k"];

    vector<vector<double>> A; // rows: [x^2, y]
    vector<double> C;         // constant term vector (all 1)

    for (auto& item : test.items()) {
        if (item.key() == "keys") continue;

        double x = stod(item.key());
        int base = stoi(item.value()["base"].get<string>());
        string val = item.value()["value"];
        double y = decodeValue(val, base);

        A.push_back({x * x, y});
        C.push_back(1); // fixed c = 1
    }

    if (k > (int)A.size()) k = (int)A.size();
    A.resize(k);
    C.resize(k);

    auto [a, b] = solveForAandB(A, C);

    cout << fixed << setprecision(6);
    cout << "a = " << a << ", b = " << b << ", c = 1\n\n";
}

int main() {
    ifstream fin("input.json");
    if (!fin) {
        cerr << "Cannot open input.json\n";
        return 1;
    }

    json j;
    fin >> j;

    if (j.is_array()) {
        for (auto& test : j) {
            processTestCase(test);
        }
    } else if (j.is_object()) {
        processTestCase(j);
    } else {
        cerr << "Input JSON must be an object or array of objects.\n";
        return 1;
    }

    return 0;
}
