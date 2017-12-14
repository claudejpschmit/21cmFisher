#pragma once

using namespace std;

class WignerPythonInterface {
    public:
        WignerPythonInterface();
        ~WignerPythonInterface();

        double W6J(int l1, int l2, int l3, int l4, int l5, int l6);
};
