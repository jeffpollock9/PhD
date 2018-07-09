#ifndef TIMER_H_
#define TIMER_H_

#include <iostream>
#include <chrono>

#define INIT_TIMER auto start = std::chrono::high_resolution_clock::now();

#define START_TIMER start = std::chrono::high_resolution_clock::now();

#define STOP_TIMER(name) std::cout << "RUNTIME of " << name << ": " << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now()-start).count() << " seconds " << std::endl;

#endif /* TIMER_H_ */
