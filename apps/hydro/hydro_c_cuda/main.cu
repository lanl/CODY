#include "hydro_args.h"
#include "engine.h"

int main(int argc,char *argv[]){
    hydro_args H;

    process_args(argc,argv,&H);

    engine(H);
}

