#include <iostream>
#include "Bigint.h"

using namespace std;

int main() {

    while (1) {
        Bigint tmp;
        Bigint tmp2;
        cin >> tmp >> tmp2;

        cout<<tmp<<"+"<<tmp2<<"="<<tmp+tmp2<<endl;
        cout<<tmp<<"-"<<tmp2<<"="<<tmp-tmp2<<endl;
        cout<<tmp<<"*"<<tmp2<<"="<<tmp*tmp2<<endl;
        cout<<tmp<<"/"<<tmp2<<"="<<tmp/tmp2<<endl;
        cout<<"pow("<<tmp<<", 2) = "<<pow(tmp, 2)<<endl;
        cout<<"pow("<<tmp<<", 3) = "<<pow(tmp, 3)<<endl;

    }

    return 0;
}
