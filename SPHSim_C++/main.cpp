//
//  main.cpp
//  SPHSim_C++
//
//  Created by Peihong Guo on 4/22/14.
//  Copyright (c) 2014 Peihong Guo. All rights reserved.
//

#include <iostream>
#include "SPHSystem.h"

int main(int argc, const char * argv[])
{
    SPHSystem *sph;
    if( argc > 1 )
        sph = new SPHSystem(argv[1]);
    else
        sph = new SPHSystem;
    
    sph->run();
    
    std::cout << "Done!\n";
    return 0;
}

