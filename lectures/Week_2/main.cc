#include<iostream>
#include<cstdio>
#include"hello.h"
#include"sfuns.h"
int main(){
 hello();
 double x=1;
 double y = sfuns::fgamma(x);
 std::cout << "fgamma(1)="<< y<< "\n";
 std::printf("fgamma(1)=%g\n",y);
 for(double x=1;x<=9;x+=1){
    std::cout <<"fgamma(" << x <<") =" << sfuns:: fgamma(x)
        << " tgamma("<< x <<")="<< std::tgamma(x) << "\n";
    }
 return 0;
}
