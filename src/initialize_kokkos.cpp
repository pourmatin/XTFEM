//
//  initialize_kokkos.cpp
//  
//
//  Created by Hossein Pourmatin on 8/22/16.
//
//

#include <stdio.h>
#include <Kokkos_Core.hpp>

extern "C" {
    
    /*void kokkos_init_(int& argc, char* argv[], int len) {
        Kokkos::initialize(argc, argv);
    }*/

    void kokkos_init_() {
        Kokkos::InitArguments args;
        // 12 (CPU) threads per NUMA region
        args.num_threads = 1;
        args.device_id = 0;
        Kokkos::initialize(args);
    }
    
    void kokkos_finish_() {
        Kokkos::finalize ();
    }
}
