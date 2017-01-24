#include "exp_funs.h"

int main(){
    for (double i=-0.5; i<=0.5; i+=0.001){
        exp_funs(i);
        //exp_funs_Taylor(i);
    }
    return 0;
}
